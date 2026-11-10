# Mixed-precision GPU numerics

## Scope of this chapter

Numerical techniques introduced by the GPU version of OpenMX (H. Kawai) to obtain FP64-quality results on GPUs whose
FP64 throughput is weak: the device precision policy, double-float grid integrals, emulated FP64 GEMM, and the
FP32 eigensolve with iterative refinement. They are written as algorithms so that they can be implemented independently
(see [GPU strategy](../design/gpu-strategy.md)). Facts observed in the GPU version are tagged [Convention];
choices for the Rust implementation are tagged [Design]. The techniques are credited to H. Kawai.
The [Convention] facts describe the GPU version at commit 08d0521 (the other LCAO chapters were checked against 3871b1f).

## Notation [Design]

| Symbol | Meaning |
|---|---|
| \\( A \\) | Dense \\( n\times n \\) real symmetric or complex Hermitian matrix in FP64 (after the Löwdin transformation) |
| \\( X \\) | \\( n\times n \\) matrix of approximate eigenvectors (columns), FP64 |
| \\( X\_1 \\) | The \\( k\_c \\) columns that are refined (\\( k\_c \\) = number of wanted states \\( m\_{\max} \\)) |
| \\( X\_c \\) | The remaining \\( n - k\_c \\) columns |
| \\( X\_{\mathcal C} \\) | The columns of one near-degenerate cluster (Section 4.4) |
| \\( e^{(0)} \\) | Eigenvalues from the FP32 solve |
| \\( f\_i \\) | Occupation of state \\( i \\) at the current chemical potential |
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

Input keywords: `scf.fp32.enable off` disables the optional FP32 and double-float paths (grid integrals, orbital
evaluation, FP32 eigensolve with refinement); it does not affect the emulated GEMM or the vendor FP64 emulation, which have
their own settings. `scf.gemmul8.enable` forces the emulated GEMM on or off (unset: the device policy).

**Inherited float products.** Upstream OpenMX forms some products in single precision and accumulates them in double
(neutral-atom projector terms, force and stress projectors, stored orbital values). The GPU version keeps this rounding
even when all optional reduced precision is off: promoting these products to FP64 shifted total energies by
\\( 1.6\times10^{-7} \\) Ha and failed the \\( 10^{-7} \\) gate against upstream. [Design] The Rust implementation
specifies this rounding explicitly per module and offers an FP64 mode for its own verification.

A plain-FP32 grid mode also exists (lossy FP32 while NormRD \\( > 10^{-5} \\), then FP64 with a restart of the mixing
history). It is off by default because it slows SCF convergence (27-40 steps instead of 23 in the test case). [Design]
Not adopted.

## 2. Double-float grid integrals [Convention]

[Public] Double-float ("double-single") arithmetic follows T. J. Dekker, Numer. Math. 18, 224-242 (1971); the
error-free transformations (TwoSum, FMA-based TwoProduct) are those of T. Ogita, S. M. Rump and S. Oishi,
SIAM J. Sci. Comput. 26, 1955-1988 (2005).

[Convention] A double-float value is an unevaluated sum \\( x = x\_{\mathrm{hi}} + x\_{\mathrm{lo}} \\) of two FP32 numbers with
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
- Workspace: in total at most 30% of device memory (50% in the noncollinear band solver), divided among the processes
  sharing the device; at least 1.5 GiB left free; and a per-process cap of 256 MiB, above which the product is computed in
  blocks rather than abandoned. If the workspace cannot be obtained, the product falls back to native FP64 GEMM.
- Products routed through it on weak-FP64 devices: the Löwdin forward and back transformations, density-matrix products,
  divide-and-conquer and Krylov products, the refinement products of Section 4, and the tile products of distributed GEMM.
- On A100 it is used only for large real forward transformations (\\( n \ge 6144 \\)), where it is 7-21% faster than
  native FP64; in a benchmark of 12 general products (real and complex, workspace capped), native FP64 was faster in 11.

Accuracy: against long-double references, the maximum relative error was \\( 3.7\times10^{-13} \\) over 39 general real
and complex products (both scaling modes). In a strict cold-start comparison (650-atom band system, `scf.criterion`
\\( 10^{-12} \\), AMD GPU), runs with the default emulation policy and with all accelerated paths in native FP64 took the
same number of SCF steps and differed by \\( 2.6\times10^{-11} \\) Ha in total energy and \\( 5\times10^{-12} \\) Ha/bohr in
forces.

## 4. FP32 eigensolve with refinement [Convention]

### 4.0 Underlying algorithm [Public]

The refinement is the iterative refinement of a symmetric eigendecomposition by T. Ogita and K. Aishima,
Japan J. Indust. Appl. Math. 35, 1007-1035 (2018); clustered eigenvalues are treated in part II of that work,
Japan J. Indust. Appl. Math. 36, 435-459 (2019). Given an approximate eigenvector matrix \\( \hat X \\) of a real
symmetric \\( A \\), one step uses \\( R = I - \hat X^{\mathsf T}\hat X \\) and \\( S = \hat X^{\mathsf T} A \hat X \\),
the Rayleigh-quotient estimates \\( \tilde\lambda\_i = s\_{ii}/(1 - r\_{ii}) \\), the threshold
\\( \delta = 2\left(\lVert S - \mathrm{diag}(\tilde\lambda)\rVert\_2 + \lVert A\rVert\_2\lVert R\rVert\_2\right) \\), and the
update \\( \hat X \leftarrow \hat X + \hat X \tilde E \\) with
\\( \tilde e\_{ij} = (s\_{ij} + \tilde\lambda\_j r\_{ij})/(\tilde\lambda\_j - \tilde\lambda\_i) \\) when
\\( \lvert\tilde\lambda\_j - \tilde\lambda\_i\rvert > \delta \\) and \\( \tilde e\_{ij} = r\_{ij}/2 \\) otherwise. For well
separated eigenvalues the error decreases quadratically; for a pair closer than \\( \delta \\) only orthonormality is
restored, and the pair's invariant subspace has to be resolved by other means. The Hermitian case replaces transposes by
conjugate transposes. Sections 4.1-4.5 record how the GPU version adapts this scheme.

### 4.1 Preparation

1. \\( A \\) is made exactly symmetric (Hermitian) by taking its lower triangle as authoritative; in the Hermitian case the
   imaginary parts of the diagonal are set to zero.
2. **Cold start**: all \\( n \\) eigenpairs of \\( A \\) rounded to FP32 are computed with a two-stage dense solver. If that
   solver does not converge, one retry with a one-stage solver is made; if that also fails, FP64 is used for the rest of the
   SCF cycle. \\( X \\) is the FP32 eigenvector matrix widened to FP64 and \\( e^{(0)} \\) the FP32 eigenvalues; any
   non-finite value rejects the solve.
3. **Warm start** (allowed from the second SCF step on, when a basis of the same order from a successful refined solve is
   stored): \\( X \\) and the eigenvalue estimates are taken from the previous solve. The warm start is accepted only if the
   quantities of the first correction step (Section 4.3) evaluated on the old basis satisfy
   \\( \max\_{i\ne j}\lvert S\_{ij}\rvert \le 2\times10^{-5} \\) and \\( \max\lvert R\_{ij}\rvert \le 2\times10^{-5} \\); otherwise a
   cold start is made. (A fresh FP32 solve leaves about \\( 4\times10^{-6} \\).)
4. Spectral radius estimate \\( a = \max(\lvert e^{(0)}\_1\rvert, \lvert e^{(0)}\_n\rvert) \\) from the current estimates.

### 4.2 Occupations

The chemical potential is fixed by bisection (200 halvings) so that the Fermi occupations of the lowest
\\( m\_{\max} \\) estimates of each spin add up to the electron count. With spin polarization the estimates of both spins
enter one common chemical potential, so the processes that own the two spins exchange their estimates first. In band
calculations the chemical potential is determined from the estimates of all k-points together. If one process holds
more than one k-point, each of its problems is refined twice: once using the chemical potential of the preceding SCF step,
and again, starting from that result, once the current chemical potential is known.

### 4.3 Correction steps

The GPU version refines only the first \\( k\_c \\) columns \\( X\_1 \\) of the full basis \\( X \\). One correction step,
applied \\( K \\) times (default \\( K = 2 \\)), is defined as follows; all matrix products use the emulated GEMM of Section 3.

- \\( S = X^\dagger A X\_1 \\) and \\( R = I\_{n\times k\_c} - X^\dagger X\_1 \\) (both \\( n\times k\_c \\)).
- Estimates \\( \lambda\_j = S\_{jj}/(1 - R\_{jj}) \\) for \\( j \le k\_c \\); the unrefined states \\( i > k\_c \\) keep their
  current estimates.
- Threshold: \\( \delta = 2\left(\max\_{i\ne j}\lvert S\_{ij}\rvert + a\max\lvert R\_{ij}\rvert\right) \\). Entrywise maxima and
  \\( a \\) replace the 2-norms of Section 4.0. \\( \delta \\) is evaluated in the first step only and kept for the remaining
  steps; after a warm start it is not allowed to be smaller than its value in the most recent cold solve.
- Correction matrix (\\( n\times k\_c \\)):
  $$
  E\_{ij} =
  \begin{cases}
  \dfrac{S\_{ij} + \lambda\_j R\_{ij}}{\lambda\_j - \lambda\_i} & \lvert\lambda\_j - \lambda\_i\rvert > \delta \text{ and } \lvert f\_i - f\_j\rvert \ge 10^{-12},\\\\[1ex]
  R\_{ij}/2 & \text{otherwise.}
  \end{cases}
  $$
  Unlike Section 4.0, pairs of equal occupation are only re-orthonormalized: the density matrix depends on the occupied
  subspace, not on the basis chosen inside it.
- **Complement rotation** (when \\( k\_c < n \\)): the remaining columns \\( X\_c \\) are corrected against \\( X\_1 \\) by the
  same rule, with \\( S\_{rc} = X\_1^\dagger A X\_c \\), \\( R\_{rc} = -X\_1^\dagger X\_c \\) and the resulting
  \\( k\_c\times(n-k\_c) \\) matrix \\( E\_{rc} \\): \\( X\_c \leftarrow X\_c + X\_1 E\_{rc} \\).
  Without this rotation, the FP32-level components of \\( X\_c \\) along the occupied subspace leave the refined columns
  converging to vectors that deviate from that subspace by about \\( 10^{-8} \\); this cost 2-5 extra SCF steps and a
  \\( 2\times10^{-10} \\) Ha bias across chains of warm starts.
- Update \\( X\_1 \leftarrow X\_1 + X E \\). Both updates of a step use the basis as it was at the start of the step.
- Decision rule: when refinement is on by default (not requested explicitly) and any of its products had to fall back to
  native FP64 GEMM, refinement is switched off for the rest of the run, since native FP64 products make it slower than an
  FP64 eigensolve.

### 4.4 Rayleigh-Ritz on near-degenerate clusters

After the correction steps all entries of \\( X\_1 \\) and \\( \lambda \\) must be finite. The estimates
\\( \lambda\_1, \dots, \lambda\_{k\_c} \\) (in index order) are partitioned into maximal runs in which neighbouring
estimates differ by at most \\( 30\delta \\). A run of two or more states whose occupations, re-evaluated at the refined
estimates, span a range of at least \\( 10^{-12} \\) is a cluster. For each cluster with columns \\( X\_{\mathcal C} \\),
the small generalized problem \\( H\_{\mathcal C} y = \mu N\_{\mathcal C} y \\) with
\\( H\_{\mathcal C} = X\_{\mathcal C}^\dagger A X\_{\mathcal C} \\) and \\( N\_{\mathcal C} = X\_{\mathcal C}^\dagger X\_{\mathcal C} \\)
is solved in FP64, and \\( X\_{\mathcal C} \\) and its estimates are replaced by the Ritz vectors and values. (The products
use the emulated GEMM when they amount to about \\( 10^9 \\) flops or more and native FP64 below; problems of order
\\( < 256 \\) are solved on the host.) If a cluster has more than \\( \min(2048, \sqrt{n k\_c/2}) \\) states, or a solve
fails, FP64 is used for the rest of the SCF cycle. On success the basis and estimates are stored for the next warm start.

### 4.5 Failure handling and SCF guards

- \\( A \\) is never overwritten before refinement succeeds, so an FP64 solve is always possible.
- Failures are classified by scope: this solve only, the rest of the SCF cycle, or the rest of the run. Failures in the
  first SCF step affect that solve only.
- The SCF step that ends a converged run is solved again in FP64 if it was refined (the energy-weighted density matrix
  from refined vectors can be off by \\( 2\times10^{-6} \\)). The last step allowed by `scf.maxIter` is always FP64.
- A stop on a refined step additionally requires NormRD \\( \le 100\times \\) `scf.criterion`. Without this guard one
  test run stopped five steps early with a \\( 2\times10^{-9} \\) Ha error.
- Refinement is not used for: a single process solving both spins (including the serialized spin mode), core-hole
  (XANES) calculations, calculations with emptied orbitals or states, or matrices below order 1000 on AMD devices.

### 4.6 Vendor-solver workaround on AMD [Convention]

The AMD dense Hermitian solver returned non-finite eigenvectors for real-valued or block-diagonal Hermitian input. The GPU
version applies a diagonal phase similarity \\( D = \mathrm{diag}(e^{2\pi i\\,\mathrm{frac}(k\varphi)}) \\) with the golden
ratio \\( \varphi \\) before the FP32 solve, and replaces exactly-zero off-diagonal entries by about \\( 10^{-25} \\).
The failure was observed with ROCm 7.2. [Open] Whether later AMD libraries still need this.

### 4.7 Measured effect

216-atom noncollinear cluster (\\( 2n = 5616 \\)) on an RTX 5080: the whole run went from 261.3 s to 197.9 s, with an
unchanged number of SCF steps; the two-stage FP32 solve took 323 ms against 638 ms for the one-stage one. On an A100 (where
refinement is not on by default), warm starts reduced the eigensolver phase from 706 to 140 ms per step. Refinement was
not adopted for the distributed (ELPA/COSMA) path: every correction step requires four or more distributed
\\( n\times n\times n \\) products, and one such product cost about as much as the FP32 solve saved.

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
