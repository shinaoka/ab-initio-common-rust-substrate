# List of public documents

These are the documents that implementers may read directly. URLs were checked as of 2026-11-05.

## Manual

- OpenMX Ver. 4.0 manual: [HTML](https://www.openmx-square.org/openmx_man4.0/index.html),
  [PDF](https://www.openmx-square.org/openmx_man4.0/openmx4.0.pdf),
  [Japanese version](https://www.openmx-square.org/openmx_man4.0jp/index.html)
- Sections particularly relevant to reimplementation:
  - 11 Functionals, 12 Basis functions (12.1 Overview, 12.2 primitive basis, 12.4 Optimization of PAO), 13 Pseudopotentials
  - 14 Cutoff energy and grid, 15 k points, 16 SCF convergence
  - 27 O(N) methods, 40 DFT+U, 45 Exchange interactions, 46 NEGF, 47 Wannier functions
  - 70 Interface for post-processing codes (70.1 scfout file), 71 OpenMX API
- ADPACK Ver. 2.2 manual: [HTML](https://www.openmx-square.org/adpack_man2.2/adpack2_2.html)

## Technical notes

Index: <https://www.openmx-square.org/tech_notes/tech_notes.html>

| Title | URL | Related chapter of this note |
|---|---|---|
| Total Energy and Forces (Ver. 1.3) | [tech1-1_3](https://www.openmx-square.org/tech_notes/tech1-1_3/tech1-1_3.html) | Two-center integrals, total energy and forces |
| Stress tensor | [tech23](https://www.openmx-square.org/tech_notes/tech23-1_0.pdf) | Total energy and forces |
| Efficient Projector Expansion in LCAO | [tech15](https://www.openmx-square.org/tech_notes/tech15-1_0/tech15-1_0.html) | Neutral-atom potential |
| Numerical Atomic Basis Orbitals from H to Kr | [tech14](https://www.openmx-square.org/tech_notes/tech14-1_0/tech14-1_0.html) | PAO basis |
| Non-Collinear Spin Density Functional | [tech2](https://www.openmx-square.org/tech_notes/tech2-1_0/tech2-1_0.html) | (not yet written) |
| O(N) Krylov Subspace Method | [tech3](https://www.openmx-square.org/tech_notes/tech3-1_1/tech3-1_1.html) | (not yet written) |
| O(N) Divide-Conquer with LNO | [tech26](https://www.openmx-square.org/tech_notes/tech26-1_0/tech26-1_0.html) | (not yet written) |
| Efficient Integration of Green's Functions | [tech4](https://www.openmx-square.org/tech_notes/tech4-1_1/tech4-1_1.html) | Density matrix |
| Orbital Optimization Method | [tech5](https://www.openmx-square.org/tech_notes/tech5-1_0/tech5-1_0.html) | PAO basis |
| LDA+U Method | [tech6](https://www.openmx-square.org/tech_notes/tech6-1_0/tech6-1_0.html) | (not yet written) |
| Closest Wannier Function | [cwf](https://www.openmx-square.org/tech_notes/cwf/cwf.html) | (not yet written) |
| Geometry Optimization | [tech7](https://www.openmx-square.org/tech_notes/tech7-1_0/tech7-1_0.html) | (not yet written) |
| Extrapolation of Charge Density | [tech8](https://www.openmx-square.org/tech_notes/tech8-1_0/tech8-1_0.html) | Density matrix and charge mixing |
| Charge Analysis | [tech9](https://www.openmx-square.org/tech_notes/tech9-1_0/tech9-1_0.html) | (not yet written) |
| Householder Method for Tridiagonalization | [tech10](https://www.openmx-square.org/tech_notes/tech10-1_0/tech10-1_0.html) | Generalized eigenvalue |
| Electric Polarization by Berry Phase | [tech11](https://www.openmx-square.org/tech_notes/tech11-1_1/tech11-1_1.html) | (not yet written) |
| Charge Mixing Methods | [tech12](https://www.openmx-square.org/tech_notes/tech12-1_0/tech12-1_0.html) | Density matrix and charge mixing |
| Two Electron Repulsion Integrals | [exx/JCP2009](https://www.openmx-square.org/exx/JCP2009/index.html) | (not yet written) |
| Fast Spherical Bessel Transform | [exx/CPC2009](https://www.openmx-square.org/exx/CPC2009/index.html) | Two-center integrals |
| Electronic Transport Calculations | [tech13](https://www.openmx-square.org/tech_notes/tech13-1_0/tech13-1_0.html) | (not yet written) |
| Eigenchannels and Current density | [tech20](https://www.openmx-square.org/tech_notes/tech20-1_0/tech20-1_0.html) | (not yet written) |
| Low-order Scaling Method | [tech16](https://www.openmx-square.org/tech_notes/tech16-1_0/tech16-1_0.html) | (not yet written) |
| Finite Element Method for Atoms | [tech17](https://www.openmx-square.org/tech_notes/tech17-1_0/tech17-1_0.html) | (not yet written) |
| Unfolding Optical Transition Weights | [PDF](https://www.openmx-square.org/tech_notes/unfolding_otw.pdf) | (not yet written) |
| Absolute Binding Energies of Core Levels | [tech24](https://www.openmx-square.org/tech_notes/tech24-1_0/tech24-1_0.html) | (not yet written) |
| Calculation of Optical Matrix Elements | [PDF](https://www.openmx-square.org/tech_notes/Opt_MatElements.pdf) | (not yet written) |
| Dielectric Function and Conductivity | [PDF](https://www.openmx-square.org/tech_notes/Dielectric_Function_YTL.pdf) | (not yet written) |
| Spin texture in k-space | [PDF](https://www.openmx-square.org/tech_notes/note_kSpin-1_0.pdf) | (not yet written) |
| Nudged Elastic Band Method | [PDF](https://www.openmx-square.org/tech_notes/NEB.pdf) | (not yet written) |

## Papers

From Section 3 "Related papers" of the manual.

| Topic | Reference |
|---|---|
| Variational optimization of PAO basis | T. Ozaki, PRB **67**, 155108 (2003) |
| Numerical atomic basis for H to Kr | T. Ozaki and H. Kino, PRB **69**, 195113 (2004) |
| Projector expansion | T. Ozaki and H. Kino, PRB **72**, 045121 (2005) |
| Fast spherical Bessel transform | J. Toyoda and T. Ozaki, CPC **181**, 277 (2010) |
| Krylov subspace O(N) method | T. Ozaki, PRB **74**, 245101 (2006) |
| Low-order scaling method | T. Ozaki, PRB **82**, 075131 (2010) |
| DC-LNO method | T. Ozaki, M. Fukuda, G. Jiang, PRB **98**, 245137 (2018) |
| Massive parallelization | T. V. T. Duy and T. Ozaki, CPC **185**, 777 (2014); CPC **185**, 153 (2014) |
| DFT+U | M. J. Han, T. Ozaki, J. Yu, PRB **74**, 045110 (2006) |
| Exchange interactions | M. J. Han, T. Ozaki, J. Yu, PRB **70**, 184421 (2004) |
| NEGF | T. Ozaki, K. Nishio, H. Kino, PRB **81**, 035116 (2010) |
| ESM method | T. Ohwaki et al., JCP **136**, 134101 (2012) |
| MLWF | H. Weng, T. Ozaki, K. Terakura, PRB **79**, 235118 (2009) |
| CWF | T. Ozaki, PRB **110**, 125115 (2024) |
| Band unfolding | C.-C. Lee, Y. Yamada-Takamura, T. Ozaki, JPCM **25**, 345501 (2013) |
| Binding energies of core levels | T. Ozaki and C.-C. Lee, PRL **118**, 026401 (2017) |
| Accuracy comparison (Δ gauge) | K. Lejaeghere et al., Science **351**, aad3000 (2016) |

For general methods (Kleinman–Bylander form, Blöchl's multiple projectors, Troullier–Martins pseudopotentials,
Morrison–Bylander–Kleinman form, Pulay/Kerker mixing, etc.), the original papers are referenced.
