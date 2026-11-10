# Licensing

## This book

- Text: [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)
- Scripts, code samples and test data: MIT License

Code samples and test data are MIT so that they can be incorporated into implementations under any license.

## Implementations

Implementations built from this book are intended to be released under MIT/Apache-2.0
(dual license, the Rust convention). This is possible only for implementations that do not copy or
adapt code from GPL or proprietary reference codes; see [Clean-room procedure](clean-room.md).

## Reference codes and data

| Reference | License | How we use it | Notes |
|---|---|---|---|
| OpenMX (code) | GPL-3.0 | Clean-room reference for the LCAO family; fixtures | See "OpenMX" below |
| OpenMX PAO/VPS databases | GPL (distributed with OpenMX) | Basis sets and pseudopotentials for the LCAO family | May be incorporated with a statement of origin; see below |
| OpenMX GPU version (H. Kawai) | GPL-3.0 | Reference for [Convention] facts and fixtures; GPU techniques specified as algorithms | Conventions checked against commit 3871b1f (the GPU chapters, against 08d0521); see [GPU strategy](../methods/lcao-openmx/design/gpu-strategy.md) for what may flow into the implementation |
| cuSOLVER call sequence in OpenMX 4.0 (also in the GPU version) | NVIDIA proprietary notice | Not used | Re-specify from public cuSOLVER documentation |
| FFTW | GPL-2.0+ | Not linked | Use permissive FFTs (rustfft, cuFFT, rocFFT) |
| ELPA | LGPL-3.0 | Distributed eigensolver via FFI | Link dynamically |
| MAGMA, COSMA, GEMMul8 | BSD-3 / BSD-3 / MIT | GPU eigensolvers, distributed GEMM, emulated FP64 GEMM via FFI | Permissive; GEMMul8 LICENSE (MIT, RIKEN R-CCS) checked; keep notices if code is reused |
| Quantum ESPRESSO, ABINIT, GPAW, Elk, exciting | GPL | Clean-room references; fixture generation | Code must not be copied; libPAW (in ABINIT) is treated as GPL |
| WIEN2k, VASP | Proprietary | Published results only | — |
| SIRIUS | BSD-3-Clause | Algorithm reference for PAW and FP-LAPW; possible verification oracle | Permissive; keep notices if code is reused |
| FLEUR | MIT | Algorithm reference and oracle for FLAPW | Permissive |
| dftatom | MIT | Radial Schrödinger/Dirac solvers | Permissive |
| ACWF verification data | CC BY 4.0 | All-electron equation-of-state references | Cite the data set |
| PseudoDojo / JTH PAW data | No license file found | — | Redistribution status unknown; do not bundle until clarified |

A fuller survey of permissive codes, libraries and data is in [Prior art](../substrate/prior-art.md).

## OpenMX

The approach was discussed with the principal author of OpenMX (T. Ozaki) in November 2026.
The points confirmed were:

- Building new software from the GPL-licensed OpenMX is not a problem.
- An **independent implementation that does not copy or adapt the existing OpenMX code** may be
  released under MIT/Apache-2.0.
- The **PAO and VPS databases** may be incorporated (they are already reused by several codes,
  including commercial ones). Their origin must be stated when they are distributed.
- To avoid confusion, implementations should use a name different from "OpenMX".
- When results are published, the related OpenMX papers should be cited:
  <https://www.openmx-square.org/openmx_man4.0jp/s3_related-papers.html>.

## Third-party documents

Technical notes, papers and manuals of reference codes have their own copyrights. This book does not
transcribe their text, figures or tables; it restates equations and facts in its own words and cites sources.

## Fixtures

Outputs of existing implementations (matrix elements, energies, ...) are numerical program outputs.
Each fixture records the program and version, the input file, and the names of the basis sets and
pseudopotentials used. Fixtures that use the OpenMX PAO/VPS databases may include those files,
together with a statement of their origin. Data whose redistribution status is unknown
(see the table above) is not stored; fixtures record how to obtain it instead.
