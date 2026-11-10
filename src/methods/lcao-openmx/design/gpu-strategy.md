# GPU strategy and porting policy [Design]

This chapter records what the GPU version of OpenMX by H. Kawai (the "GPU version"; analysed at its November 2026 state)
teaches us, and how its results are carried into the Rust implementation. It describes techniques and policy only.
Following the [writing rules](../../../process/writing-rules.md), it names no source files, functions or variables of
the GPU version; numerical techniques are specified in [Mixed-precision GPU numerics](../modules/gpu-mixed-precision.md).

## 1. What the GPU version accelerates [Convention]

The GPU version keeps the program structure of OpenMX and offloads selected stages. Most offload is done with
compiler directives (OpenACC for NVIDIA, OpenMP target offload for AMD) inside the existing C routines; a small number
of hand-written kernels cover the grid integrals, and dense linear algebra goes to vendor and third-party libraries.

| Stage | Where it runs | Technique |
|---|---|---|
| Orbital values on the real-space grid | GPU (up to \\( l \le 3 \\) and 8 radial functions per \\( l \\)), else host SIMD | Directive loops; optionally evaluated on the fly inside Hamiltonian and density tiles |
| Overlap of atom spheres with grid cells | GPU, else host SIMD | Directive loops |
| Electron density on the grid from the density matrix | GPU | Hand-written kernel, one output point per thread, sparse pair lists |
| Grid matrix elements \\( \langle\phi\lvert \delta V\_H + V\_{xc} + V\_{NA}\rvert\phi\rangle \\) | GPU | Hand-written tiled kernel: one thread block per atom pair, grid points staged in shared memory |
| Overlap and kinetic two-center integrals | CPU (GPU is transfer-bound) | Cached spherical Bessel transforms, host SIMD |
| Nonlocal projector overlaps, neutral-atom projector expansion | GPU | Directive loops; Bessel tables built on the device, Gaunt coefficients and spherical harmonics on the host |
| XC on the grid (GGA) | GPU | Gradient stencil, functional, divergence, noncollinear rotation |
| Poisson equation (FFT) | CPU | Host FFT library |
| Grid sums of the total energy | GPU | Directive reductions |
| Dense eigenproblems (Cluster, Band, divide-and-conquer, DC-LNO, Krylov) | GPU above a size threshold, else CPU | Vendor eigensolvers; Löwdin transformation by GEMM; optional FP32 solve with refinement |
| Density matrix and energy-weighted density matrix | GPU (eigenvectors stay on the device) | Directive and hand-written reductions |
| Charge mixing (RMM-DIISH) | GPU Gram matrix, CPU otherwise | Directive reduction plus one collective |
| Forces | Mostly GPU | Directive loops |
| Stress, NEGF | CPU | — |

Size thresholds below which dense problems stay on the CPU: collinear cluster \\( n < 400 \\), noncollinear cluster
\\( 2n < 400 \\), band \\( n < 800 \\) (AMD: 400), divide-and-conquer 400 (AMD: 300), DC-LNO 800 (AMD: 400), Krylov 300.

### Measured effect [Convention]

Against upstream OpenMX on the CPU of the same machine:

| Hardware | 216-atom noncollinear cluster | 650-atom band system | Upstream "L" test suite |
|---|---|---|---|
| RTX 5080 (18 MPI ranks sharing the GPU) | 8.2× | 5.7× | 1.7× |
| A100 40 GB | 3.9× | 5.0× | 1.3× |
| H100 PCIe 80 GB | 4.4× | 4.5× | 1.4× |
| MI300A (4 APUs) | 1.9× | 2.5× | 0.7× |

Small systems are slower on every GPU (launch latency and grid-integral overhead). On data-centre GPUs dense
diagonalization remains the dominant phase; on consumer GPUs the gains come largely from reduced and emulated precision.

### Validation practice [Convention]

The GPU version compares total energies and **every force component** (never force sums), together with SCF iteration
counts, pairwise between upstream on the CPU, the GPU version built for the CPU, and the GPU version on the GPU.
The pass gate is \\( 10^{-7} \\) Ha and \\( 10^{-7} \\) Ha/bohr. Strict tests at an SCF criterion of \\( 10^{-12} \\)
agree to about \\( 10^{-10} \\) Ha. We adopt the same gates and the same pairwise scheme
(see [Testing strategy](../../../testing/strategy.md)), with the GPU version as an additional oracle.

## 2. Provenance and what may flow into the Rust code [Design]

Almost all of the GPU version is derived from upstream OpenMX (GPL-3.0):

- Directive-based offload lives inside upstream routines.
- The C++ core library of the GPU version reorganizes upstream numerical routines, deliberately to the bit,
  rather than redesigning them.
- The hand-written kernels are new code but compute upstream's integrals on upstream's data layouts.

Consequences:

| Part of the GPU version | Route into the Rust implementation |
|---|---|
| Physics and numerics inherited from upstream (including their GPU forms) | Specification chapters of this book, written under the clean-room procedure; implemented anew |
| Numerical techniques introduced by the GPU version (mixed precision, eigenvector refinement, GEMM emulation policy, capacity fallbacks) | Specified as algorithms in [Mixed-precision GPU numerics](../modules/gpu-mixed-precision.md), credited to H. Kawai; implemented anew |
| Thin bridges to external libraries, build scripts, library patches written by H. Kawai | May be reused directly if H. Kawai agrees to release them under a permissive license; otherwise re-specified |
| Files carrying an NVIDIA proprietary notice (a cuSOLVER call sequence inherited from upstream 4.0) | Not used; the call sequence is re-specified from public cuSOLVER documentation |
| Test fixtures that embed upstream code verbatim | Not used |

Algorithms are not protected by copyright; their expression in code is. Writing the techniques down as algorithms keeps the
Rust code independent while preserving the know-how. Because H. Kawai has read the upstream source extensively, his natural
roles in the clean-room procedure are spec writer and reviewer for GPU techniques, and owner of the validation oracle.

## 3. External libraries [Design]

All of the following expose C (or ScaLAPACK-style) entry points and can be called from Rust through FFI.
They are wrapped behind traits in the substrate so that each method can switch implementations.

| Library | License | Used for | FFI notes |
|---|---|---|---|
| cuBLAS, cuSOLVER (NVIDIA) | NVIDIA EULA (runtime linking) | GEMM; dense symmetric/Hermitian eigensolvers, including the emulated-FP64 modes of recent cuSOLVER | Device pointers, column-major. Accessed through tenferro-rs |
| hipBLAS, hipSOLVER, rocSOLVER (AMD) | MIT | Same on AMD | Same |
| GEMMul8 (RIKEN R-CCS) | MIT (to be confirmed from its LICENSE file) | FP64 GEMM emulated with INT8 tensor cores (Ozaki scheme II) | C++ templates: needs a small C shim compiled with nvcc/hipcc; link as a prebuilt library. Use public entry points only |
| MAGMA | BSD-3 | Single- and multi-GPU dense eigensolvers | C API |
| ELPA | LGPL-3.0 | Distributed eigensolvers on CPU and GPU | C API; link dynamically. Expects a Fortran MPI communicator handle and BLACS block-cyclic distribution |
| COSMA (+ COSTA, Tiled-MM) | BSD-3 | Communication-optimal distributed GEMM | ScaLAPACK-compatible C entry points; avoid internal C++ APIs |
| ScaLAPACK / BLACS | BSD | Block-cyclic distribution for ELPA and COSMA | Avoid non-public BLACS functions |
| FFTW | GPL-2.0+ | — | Not linked. Use a permissive FFT (rustfft on the host, cuFFT/rocFFT via tenferro-rs) |
| Kokkos | Apache-2.0 | — | Not needed in Rust |

Library patches made in the GPU version (device mapping in MAGMA, a GEMM hook in Tiled-MM, routing MAGMA's BLAS-3 phases
through GEMMul8 on AMD) are re-specified as requirements rather than carried over, unless released permissively.

## 4. Kernels in the Rust implementation [Design]

- Host code is Rust. Device kernels are CUDA C source compiled at run time with NVRTC through tenferro-rs
  (hipRTC for AMD, sharing the kernel source where possible).
- The directive loops of the GPU version have no counterpart in Rust: each becomes an explicit kernel written from the
  specification of the corresponding module.
- Kernel templates are limited to small parameter sets (real/complex, precision mode), which NVRTC handles well.
- Device data structures are flat (pair-offset tables, compressed sparse rows, column-major dense matrices, interleaved
  complex). The host-side packing from nested per-atom arrays is redesigned around the
  [sparse block format](sparse-block-format.md).
- Reproducibility contracts that matter for agreement with the oracle must be stated explicitly in kernel specifications:
  summation order where the oracle fixes it, the point where float products are widened to double, and FMA contraction
  (NVRTC options).

## 5. Process and memory model [Design]

The GPU version runs many MPI ranks per GPU (sharing through NVIDIA MPS) and therefore needs extensive machinery for memory
negotiation: grouping ranks by device, collective capacity verdicts, a registry of phase memory needs, waits for peers,
and fallback ladders. The Rust implementation uses **one process (or one host thread pool) per GPU**, with MPI across
GPUs and nodes. Most of the negotiation machinery then disappears; what remains is specified as:

- a per-phase memory budget computed from problem sizes and library workspace queries;
- a capacity ladder for dense eigenproblems: single GPU → multi-GPU on one node → distributed (ELPA/COSMA) → CPU;
- a size threshold below which a stage stays on the CPU.

## 6. Phases [Design]

| Phase | Content | Done when |
|---|---|---|
| 0 | CPU reference implementation in Rust from the specifications | Small-system fixtures pass at the [Convention] tolerances; agrees with the GPU version's CPU build within \\( 10^{-7} \\) |
| 1 | Dense linear algebra on GPU through vendor libraries (Löwdin transformation, eigensolvers, density matrix with eigenvectors kept on the device) | Same gates on the GPU; diagonalization speedup measured on a 216-atom noncollinear cluster |
| 2 | Grid kernels via NVRTC: orbital evaluation, density, grid matrix elements, XC | Same gates; per-kernel timings compared with the GPU version |
| 3 | Mixed precision: GEMMul8 through FFI, FP32 eigensolve with refinement, double-float grid integrals | Strict test (\\( 10^{-10} \\) Ha) on consumer GPUs |
| 4 | Scale-out: multi-GPU dense solvers (MAGMA), distributed solvers (ELPA, COSMA), AMD via hipRTC | Benchmarks on A100/H100/MI300A comparable to the GPU version |

## [Open]

- Whether H. Kawai will release the library bridges and patches under a permissive license.
- GEMMul8's license file (only its README has been checked).
- The tenferro-rs plan for hipRTC and for exposing cuSOLVER emulation modes.
- Whether to keep the size thresholds of the GPU version or re-tune them for the Rust implementation.
