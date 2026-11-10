# Prior art

Codes, libraries and data that inform the substrate. Licenses were checked against primary sources
(LICENSE files, repository APIs, official sites) in November 2026; re-check before reusing any code.

## Codes with permissive licenses

| Name | License | Language | What it implements | GPU | Use for us |
|---|---|---|---|---|---|
| [SIRIUS](https://github.com/electronic-structure/SIRIUS) | BSD-3-Clause | C++ | NC/USPP/PAW plane waves and full-potential LAPW in one library; backend of QE and CP2K | CUDA, ROCm | Primary design and algorithm reference for PAW and FLAPW; component-level verification oracle |
| [FLEUR](https://github.com/JuDFTteam/FLEUR) | MIT | Fortran | FLAPW with local orbitals; ACWF participant | not verified | Reference implementation and oracle for FLAPW |
| [JuKKR](https://github.com/JuDFTteam/JuKKR) | MIT | Fortran | Full-potential all-electron KKR | not verified | All-electron comparison; radial techniques |
| [DFTK.jl](https://github.com/JuliaMolSim/DFTK.jl) | MIT | Julia | Plane waves, norm-conserving; automatic differentiation | yes | SCF, mixing and differentiability design for plane waves (no PAW found) |
| [PWDFT-DCU](https://github.com/USTC-Hu-Group/PWDFT-DCU) | BSD-3-Clause | C++ | Plane waves (NC), LOBPCG/PPCG, ACE, ISDF | HIP | Eigensolvers and hybrid functionals |
| [NWChem](https://github.com/nwchemgit/nwchem) (`src/nwpw/paw`) | ECL-2.0 (Apache-2.0 derivative) | Fortran | Plane-wave PAW | partial | A readable PAW implementation under a permissive license; check ECL-2.0 terms before reuse |
| [FlapwMBPT](https://github.com/FermiQ/FlapwMBPT) | BSD-3-like (BNL) | Fortran | FLAPW + GW | not verified | LAPW reference; contains an "unless otherwise stated" clause, check per file |

## Building blocks

| Name | License | Use for us |
|---|---|---|
| [dftatom](https://github.com/certik/dftatom) | MIT | Radial Schrödinger, Dirac and Poisson solvers and atomic SCF; prime candidate for a Rust port |
| [SpFFT](https://github.com/eth-cscs/SpFFT), [SpLA](https://github.com/eth-cscs/spla), [COSTA](https://github.com/eth-cscs/COSTA) | BSD-3-Clause | Sparse plane-wave FFT and distributed matrices (GPU) |
| [spglib](https://github.com/spglib/spglib), [moyo](https://github.com/spglib/moyo) | BSD-3 / MIT OR Apache-2.0 | Symmetry; moyo is pure Rust |
| [ELSI](https://gitlab.com/elsi_project/elsi_interface) | BSD-3-Clause | Unified eigensolver interface (bundled ELPA is LGPL) |
| [ChASE](https://github.com/ChASE-library/ChASE), [MAGMA](https://github.com/icl-utk-edu/magma) | BSD-3-Clause | Large-scale (GPU) eigensolvers |
| [libgridxc](https://gitlab.com/siesta-project/libraries/libgridxc) | BSD-3-Clause | XC on radial and real-space grids |
| [sphericart](https://github.com/lab-cosmo/sphericart) | MIT | Real spherical harmonics (GPU) |
| [libcint](https://github.com/sunqm/libcint) | Apache-2.0 | Gaussian integrals (relevant to localized bases) |
| Rust crates: `rustfft`, `ndrustfft`, `lebedev_laikov`, `wigner-symbols`, `wigners` | MIT / Apache-2.0 | FFT, angular quadrature, Clebsch-Gordan and 3j symbols (Gaunt coefficients) |

## Weak copyleft and other licenses

- **libxc** (MPL-2.0): the de facto standard XC library; linking it unmodified is unproblematic. Also XCFun, libpspio (MPL-2.0).
- **INQ** (MPL-2.0, GPU-centric C++): bundles a libPAW interface, so GPL code may be mixed in.
- **libPAW** (ABINIT): its COPYING file says MPL-2.0 but source headers say GPL. **Treat as GPL**: read for understanding only, do not port.
- **ABACUS** (LGPL-3.0): plane waves and numerical atomic orbitals; a reference for OpenMX-like localized bases.
- **DFT-FE** (LGPL-2.1), **ELPA**, **DFTB+**, **simple-dftd3** (LGPL-3), **Wannier90** (LGPL-2.1).

## GPL or proprietary (papers only, or clean room)

Quantum ESPRESSO, ABINIT, GPAW, Elk, exciting, SIESTA, CP2K, SPARC, Octopus, atompaw, pslibrary,
GBRV (GPL); ecalj (AGPL-3, with a clause extending AGPL to AI-generated output); REST (GPL-2.0, Rust);
VASP and WIEN2k (proprietary). Their published methods are fair game; their code is not.

Methodologically notable among them:

- **CP2K (GPW/GAPW):** localized Gaussian basis for orbitals with plane waves for the density;
  GAPW adds PAW-like atomic augmentation to a localized basis. A precedent for LCAO + augmentation.
- **GPAW:** PAW with plane-wave, real-space-grid and LCAO bases in one code.

## Verification data

| Data | License | Use |
|---|---|---|
| [ACWF verification data](https://archive.materialscloud.org/record/2023.81) | CC BY 4.0 | All-electron (WIEN2k/FLEUR) equation-of-state references; the main external oracle |
| acwf-verification-scripts, aiida-common-workflows | MIT | Recomputation and comparison scripts |
| [SSSP 1.3](https://archive.materialscloud.org/record/2023.65) | CC BY 4.0 | Plane-wave pseudopotential library; whether each potential keeps its original license is not stated |
| pymatgen, masci-tools | MIT | I/O and analysis (masci-tools covers FLEUR I/O) |

## Not verified

- Redistribution terms of PseudoDojo/JTH PAW data, GPAW setups and the PAW-XML specification data (no license found).
- Questaal (no license statement; access on request), SG15, Delta-test data.
- GPU support of FLEUR and JuKKR.
- Files excepted by "unless otherwise stated" in FlapwMBPT.
