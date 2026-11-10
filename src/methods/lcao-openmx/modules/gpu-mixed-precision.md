# Mixed-precision GPU numerics

## Scope of this chapter

Numerical techniques introduced by the GPU version of OpenMX (H. Kawai) to obtain FP64-quality results on GPUs whose
FP64 throughput is weak: the device precision policy, double-float grid integrals, emulated FP64 GEMM, and the
FP32 eigensolve with iterative refinement. They are written as algorithms so that they can be implemented independently
(see [GPU strategy](../design/gpu-strategy.md)). Facts observed in the GPU version are tagged [Convention];
choices for the Rust implementation are tagged [Design]. The techniques are credited to H. Kawai.

## Notation [Design]

| Symbol | Meaning |
|---|---|
| \\( A \\) | Dense \\( n\times n \\) real symmetric or complex Hermitian matrix in FP64 (after the Löwdin transformation) |
| \\( X \\) | \\( n\times n \\) matrix of approximate eigenvectors (columns), FP64 |
| \\( X\_1 \\) | The \\( k\_c \\) columns that are refined (\\( k\_c \\) = number of wanted states \\( m\_{\max} \\)) |
| \\( X\_c \\) | The remaining \\( n - k\_c \\) columns |
| \\( e^{(0)} \\) | Eigenvalues from the FP32 solve |
| \\( \rho\_{\mathrm{dev}} \\) | Ratio of FP32 to FP64 throughput of the device |
| NormRD | Residual norm of charge mixing ([Density matrix and charge mixing](density-mixing.md)) |

## 1. Device precision policy [Convention]

All reduced-precision paths are keyed to one number, \\( \rho\_{\mathrm{dev}} \\). The device counts as **weak in FP64**
when \\( \rho\_{\mathrm{dev}} \ge 8 \\) (consumer GPUs such as GeForce and Radeon, workstation RTX PRO); data-centre GPUs
(A100, H100, MI300A) have \\( \rho\_{\mathrm{dev}} \le 2 \\). On NVIDIA devices the ratio is a device attribute; on AMD it
comes from a table by architecture family (data-centre CDNA parts 2, consumer parts 16).

| Place | Technique | Default |
|---|---|---|
| Grid integrals (matrix elements, density) | Double-float arithmetic (Section 2) | On when weak in FP64; FP64 otherwise |
| On-the-fly orbital evaluation | Double-float tables; the result is rounded to the FP32 value nearest the FP64 value | Same |
| Dense GEMM | Emulated FP64 by INT8 (Section 3) | On when weak in FP64; native FP64 otherwise (one exception noted in Section 3) |
| Dense eigensolver | FP32 solve with refinement (Section 4) | On when weak in FP64 **and** the emulated GEMM of the same kind (real or complex) is on |
| FP64 eigensolvers from the vendor library | Fixed-point emulation of FP64 where available, retried in native FP64 on failure | On for weak-FP64 NVIDIA devices |

Input keywords: `scf.fp32.enable off` disables every optional reduced-precision path; `scf.gemmul8.enable` forces the
emulated GEMM on or off (unset: the device policy).

**Inherited float products.** Upstream OpenMX forms some products in single precision and accumulates them in double
(neutral-atom projector terms, force and stress projectors, stored orbital values). The GPU version keeps this rounding
even when all optional reduced precision is off: promoting these products to FP64 shifted total energies by
\\( 1.6\times10^{-7} \\) Ha and failed the \\( 10^{-7} \\) gate against upstream. [Design] The Rust implementation
specifies this rounding explicitly per module and offers an FP64 mode for its own verification.

A plain-FP32 grid mode also exists (lossy FP32 while NormRD \\( > 10^{-5} \\), then FP64 with a restart of the mixing
history). It is off by default because it slows SCF convergence (27-40 steps instead of 23 in the test case). [Design]
Not adopted.

## 2. Double-float grid integrals [Convention]

A double-float value is an unevaluated sum \\( x = x\_{\mathrm{hi}} + x\_{\mathrm{lo}} \\) of two FP32 numbers with
\\( \lvert x\_{\mathrm{lo}}\rvert \le \tfrac12 \mathrm{ulp}(x\_{\mathrm{hi}}) \\). Products use an FMA-based error-free
transformation, sums use TwoSum; one term of a grid integral costs about 10 FP32 operations. Within a tile of grid points,
partial sums are kept in double-float; sums across tiles and atomic accumulations are in FP64.

Measured on an RTX 5080: matrix elements 1.42 → 0.61 s, orbital evaluation 0.77 → 0.41 s, density 0.37 → 0.098 s per rank
and SCF step, with total energies within \\( 10^{-12} \\) Ha of the FP64 kernels.

[Design] Kernel specifications must state the FMA contraction setting of the compiler (NVRTC options), since the
error-free transformations depend on it.

## 3. Emulated FP64 GEMM [Convention]

The emulated GEMM uses the Ozaki scheme II as implemented in GEMMul8 (INT8 tensor cores with the Chinese remainder theorem).

- Real products use the double-precision variant; complex products use three real INT8 products per modulus
  (about three times the cost and workspace).
- Number of moduli: 15 by default (range 2-20), accurate scaling by default.
- Workspace: at most 30% of device memory in total, at least 1.5 GiB left free, at most 256 MiB per process sharing a
  device; above the cap the product is blocked rather than abandoned. If the workspace cannot be obtained, the product
  falls back to native FP64 GEMM.
- Products routed through it on weak-FP64 devices: the Löwdin forward and back transformations, density-matrix products,
  divide-and-conquer and Krylov products, the refinement products of Section 4, and the tile products of distributed GEMM.
- On A100 it is used only for large real forward transformations (\\( n \ge 6144 \\)), where it is 7-21% faster than
  native FP64; native FP64 is faster for 11 of 12 other product types.

Accuracy checks against long-double references: maximum relative error \\( 3.7\times10^{-13} \\) over 39 general real and
complex cases; whole-program regressions changed total energies by at most \\( 1.8\times10^{-10} \\) Ha and forces by
\\( 9.6\times10^{-9} \\) Ha/bohr, with identical SCF iteration counts.

## 4. FP32 eigensolve with refinement [Convention]

### 4.1 Preparation

1. Symmetrize \\( A \\) from its lower triangle (Hermitian case: also zero the imaginary part of the diagonal).
2. **Warm start** (allowed from the second SCF step on, when a valid basis of the same order is stored): compute the
   first-iteration products \\( S = X^\dagger A X\_1 \\) and \\( G = X^\dagger X\_1 \\). If
   \\( \max\_{i\ne j}\lvert S\_{ij}\rvert > 2\times10^{-5} \\) or \\( \max\lvert I - G\rvert > 2\times10^{-5} \\), discard them
   and take the cold path. (A fresh FP32 solve leaves about \\( 4\times10^{-6} \\).)
3. **Cold start**: convert \\( A \\) to FP32 and compute all \\( n \\) eigenpairs with a two-stage dense solver; on failure
   retry once with a one-stage solver; on a second failure use FP64 for the rest of the SCF cycle. Set \\( X \\) to the FP32
   eigenvectors widened to FP64 and \\( e^{(0)} \\) to the FP32 eigenvalues. Reject any non-finite value.
4. Spectral radius \\( a = \max(\lvert e^{(0)}\_1\rvert, \lvert e^{(0)}\_n\rvert) \\).

### 4.2 Occupations

Find the chemical potential by bisection (200 iterations) so that the occupations \\( f(e^{(0)}\_i) \\) of the lowest
\\( m\_{\max} \\) states sum to the electron count. With spin polarization, the eigenvalue estimates of both spins are
exchanged first so that the chemical potential is common. In band calculations the estimates of all k-points are gathered
before the chemical potential is set; when a process holds several k-points, a first pass refines with the previous step's
chemical potential and a second pass refines again with the new one.

### 4.3 Correction steps

Repeat \\( K \\) times (default \\( K = 2 \\)):

1. \\( Y = A X\_1,\ S = X^\dagger Y,\ G = X^\dagger X\_1 \\), all with the emulated GEMM (Section 3). If an emulated product
   silently fell back to native GEMM, abandon refinement for the rest of the run.
2. \\( \lambda\_j = S\_{jj}/G\_{jj} \\) for the refined columns; \\( R = I - G \\).
3. On the first step only, set \\( \delta = 2\left(\max\_{i\ne j}\lvert S\_{ij}\rvert + a\max\lvert R\rvert\right) \\) and keep
   it fixed for the remaining steps. After a warm start, \\( \delta \\) is not allowed to fall below the value of the last cold
   solve.
4. Correction matrix:
   $$
   E\_{ij} =
   \begin{cases}
   \dfrac{S\_{ij} + \lambda\_j R\_{ij}}{\lambda\_j - \lambda\_i} & \lvert\lambda\_j - \lambda\_i\rvert > \delta \text{ and } \lvert f\_i - f\_j\rvert \ge 10^{-12},\\\\[1ex]
   R\_{ij}/2 & \text{otherwise.}
   \end{cases}
   $$
   The second case only re-orthonormalizes; rotations inside blocks of equal occupation are deliberately not resolved.
5. **Complement rotation** (when \\( k\_c < n \\)): with \\( S\_{rc} = Y^\dagger X\_c \\) and \\( G\_{rc} = X\_1^\dagger X\_c \\),
   build \\( E\_{rc} \\) by the same formula and update \\( X\_c \leftarrow X\_c + X\_1 E\_{rc} \\). Without it, the refined
   columns converge to a fixed point tilted out of the occupied subspace by about \\( 10^{-8} \\), which costs 2-5 extra SCF
   steps and a \\( 2\times10^{-10} \\) Ha bias across warm starts.
6. \\( X\_1 \leftarrow X\_1 + X E \\).

### 4.4 Rayleigh-Ritz on near-degenerate clusters

After the correction steps, check that all entries are finite. Walk the sorted \\( \lambda \\) and chain consecutive values
closer than \\( 30\delta \\) into clusters. For each cluster of size \\( \ge 2 \\) that contains different occupations,
form \\( H\_c = X\_c^\dagger A X\_c \\) and \\( N\_c = X\_c^\dagger X\_c \\) (emulated GEMM above \\( 10^9 \\) flops, native
FP64 below), solve the generalized eigenproblem (on the host for order \\( < 256 \\), on the device otherwise) and rotate
the cluster onto its solution. If a cluster exceeds \\( \min(2048, \sqrt{n k\_c/2}) \\) states, or the solve fails, use FP64
for the rest of the SCF cycle. On success, store the basis and eigenvalues for the next warm start.

### 4.5 Failure handling and SCF guards

- \\( A \\) is never overwritten before refinement succeeds, so an FP64 solve is always possible.
- Failures are classified by scope: this solve only, the rest of the SCF cycle, or the rest of the run. Failures in the
  first SCF step affect that solve only.
- The SCF step that ends a converged run is solved again in FP64 if it was refined (the energy-weighted density matrix
  from refined vectors can be off by \\( 2\times10^{-6} \\)). The last step allowed by `scf.maxIter` is always FP64.
- A stop on a refined step additionally requires NormRD \\( \le 100\times \\) `scf.criterion`. Without this guard one
  test run stopped five steps early with a \\( 2\times10^{-9} \\) Ha error.
- Refinement is not used for: a single process solving both spins, core-hole (XANES) calculations, calculations with
  emptied orbitals or states, or matrices below order 1000 on AMD devices.

### 4.6 Vendor-solver workaround on AMD [Convention]

The AMD dense Hermitian solver returned non-finite eigenvectors for real-valued or block-diagonal Hermitian input. The GPU
version applies a diagonal phase similarity \\( D = \mathrm{diag}(e^{2\pi i\\,\mathrm{frac}(k\varphi)}) \\) with the golden
ratio \\( \varphi \\) before the FP32 solve, and replaces exactly-zero off-diagonal entries by about \\( 10^{-25} \\).
[Open] Whether current AMD libraries still need this.

### 4.7 Measured effect

On an RTX 5080, the 216-atom noncollinear cluster (\\( 2n = 5616 \\)) took 197.9 s instead of 261.3 s with the same SCF step
count. On an A100, warm starts reduced the eigensolver phase from 706 to 140 ms per step; the two-stage FP32 solve took
323 ms instead of 638 ms. In distributed solves, data redistribution dominates and refinement was not adopted.

## Requirements for implementations [Design]

- The device policy (Section 1) is a single function of device properties and input keywords, evaluated once and logged.
- Every reduced-precision path has an FP64 counterpart selectable at run time, used by the verification suite.
- The refinement is implemented as a component with the inputs \\( (A, X, e^{(0)}, f, K) \\) and outputs
  \\( (X, \lambda, \text{status}) \\), independent of how \\( A \\) was produced, so that it can serve other method families.
- The emulated GEMM is reached through the substrate's GEMM trait; refinement must detect fallbacks to native GEMM.

## Tests [Design]

| Layer | Test | Tolerance |
|---|---|---|
| Analytic | Double-float sum and product against exact rational results | exact to \\( 2^{-48} \\) relative |
| Invariant | Refined vectors: residual \\( \lVert A x - \lambda x\rVert \\) and \\( \lVert X^\dagger X - I\rVert \\) by occupation class | \\( 10^{-10} \\) |
| Invariant | Random symmetric matrices with planted clusters of near-degenerate eigenvalues: refined occupied subspace equals the FP64 one | \\( 10^{-10} \\) (principal angles) |
| Comparison | Emulated GEMM against long-double reference | \\( 10^{-12} \\) relative |
| Regression | Total energy and every force component, reduced-precision mode against FP64 mode and against the GPU version | \\( 10^{-7} \\) Ha, \\( 10^{-7} \\) Ha/bohr; strict SCF \\( 10^{-10} \\) Ha |
| Fault injection | Forced solver failure, non-finite values, large clusters, emulated-GEMM fallback: correct FP64 fallback scope | — |

## [Open]

- Whether to adopt the AMD phase-similarity workaround (Section 4.6).
- Whether double-float should also be used for the force grid terms.
- Tuning of \\( K \\), the warm-start threshold and the cluster size limit for the Rust implementation.
