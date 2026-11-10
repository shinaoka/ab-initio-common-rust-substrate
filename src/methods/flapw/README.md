# FLAPW family (all-electron) [Open]

Scoping only; no specification chapters yet. This family is planned last because its details
(linearization energies, local orbitals, force corrections) are the heaviest.

## Purpose

- An all-electron reference **inside the same framework**, so that LCAO and PAW results can be verified
  against all-electron results computed with the same structure input and observables output.
- Accuracy-critical properties: electric field gradients, hyperfine fields, core-level spectroscopy, magnetism.
- GPU-native FLAPW is still rare.

## Relation to PAW

PAW and FLAPW both split space into atomic spheres and an interstitial region:

| | PAW | FLAPW |
|---|---|---|
| Inside spheres | On-site quantities expanded in radial functions × \\( Y\_{lm} \\) | Full-potential radial × \\( Y\_{lm} \\) expansion |
| Outside spheres | Plane waves (pseudo wavefunctions) | Plane waves |
| Matching | Compensation charges (multipole matching) | Value and slope matching; Weinert's pseudo-charge method |

Shared leaf components therefore include radial solvers (Schrödinger, scalar-relativistic, Dirac),
products of \\( Y\_{lm} \\)-expanded quantities via Gaunt coefficients, XC on Lebedev grids, and multipole-based
Poisson solvers. Unlike plane-wave PAW, FLAPW builds explicit dense H and S matrices, so its eigensolver usage
resembles the LCAO family.

## References to study

- SIRIUS (BSD-3-Clause): FP-LAPW alongside PAW in one GPU library.
- FLEUR (MIT): reference implementation; source of ACWF all-electron references.
- dftatom (MIT): radial solvers.
- WIEN2k (proprietary), Elk, exciting (GPL): published methods only.

See [Prior art](../../substrate/prior-art.md) for licenses.
