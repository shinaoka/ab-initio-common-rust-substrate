# Gaps in public documents

A list of items for which public documents (technical notes, papers, manual) alone are insufficient for reimplementation.
For each item, the spec writer supplements it as [Convention], or this note redefines it as [Design].

## Conventions affecting the whole

| Item | Status | Resolution |
|---|---|---|
| Ordering and signs of real spherical harmonics | The manual only says "real spherical harmonics" | Documented as [Convention] in [PAO basis](../modules/pao-basis.md) |
| Format of PAO and VPS files | Only an overview in the manual | Documented in [PAO basis](../modules/pao-basis.md) and [nonlocal pseudopotential](../modules/nonlocal-pseudopotential.md) |
| Wave-number grid and number of integration points for two-center integrals | Only examples in technical note 15 | Documented in [two-center integrals](../modules/two-center-overlap-kinetic.md) |
| Index conventions for orbitals, atoms, periodic images, and spins | Not documented | [Design] in [sparse block format](../design/sparse-block-format.md); the conventions of the existing implementation are also documented |
| How the real-space grid is determined (number of grid points from the cutoff) | Overview in technical note 1 (with an error in a coefficient) | Documented in [grid](../modules/grid-poisson-xc.md) |
| Generation of k points (including TRIM points) | Only Section 15 of the manual | Documented in [k points and Bloch sums](../modules/kpoints-bloch.md) |

## Numerical methods

| Item | Status |
|---|---|
| Thresholds and schemes for regularization of S (eigenvalue replacement, dimension reduction) | A per-path correspondence table is documented in [generalized eigenvalue](../modules/generalized-eigen.md) |
| auto-Kerker, automatic SCF tuning (manual 16.2, 16.3) | Documented in [density matrix and charge mixing](../modules/density-mixing.md) |
| Correction of the XC energy by the atomic density (enabled by default) | Not in public documents. Documented in [grid](../modules/grid-poisson-xc.md) |
| EF, RF, DIIS for geometry optimization (manual 18) | Documented in [geometry optimization and MD](../modules/geometry-md.md) (defaults differ from the manual) |
| Thermostats and barostats for MD (NVT_NH, NPT, etc.) | Documented in [geometry optimization and MD](../modules/geometry-md.md) |
| How EXX is incorporated into SCF | Unreachable in the existing implementation. Documented in [EXX and dispersion corrections](../modules/exx-vdw.md) |
| DFT-D2 and D3 | Documented in [EXX and dispersion corrections](../modules/exx-vdw.md) |
| Local problems and chemical potential of O(N) methods | Documented in [O(N) methods](../modules/on-methods.md) |

## Physical features

| Item | Status |
|---|---|
| Treatment of fully relativistic, j-dependent pseudopotentials | Documented in [nonlocal pseudopotential](../modules/nonlocal-pseudopotential.md) |
| Energy, forces, and stress with noncollinear spin and spin-orbit coupling | Documented in [noncollinear](../modules/noncollinear-energy-force.md) (stress is unsupported even in the existing implementation) |
| Constrained DFT, MAE, Zeeman term, orbital magnetic moment | Documented in [DFT+U and magnetism](../modules/dftu-magnetism.md) |
| NEGF (including noncollinear) | Documented in [NEGF](../modules/negf.md) |
| Electric field, charge doping, virtual atoms, core levels | Documented in [ESM, electric field, and core levels](../modules/esm-field-corehole.md) |
| Paths added in the GPU version | Those that do not change results are out of scope ([parallelization](../design/parallelization.md)) |

## Output and interfaces

| Item | Status |
|---|---|
| Format of the scfout file (manual 70.1) | Documented in [analysis](../modules/analysis.md) |
| OpenMX API (manual 71) | Documented in [analysis](../modules/analysis.md) |
| Mulliken, Voronoi, and ESP charges, STM, work function, etc. | Documented in [analysis](../modules/analysis.md) |

## Not yet started

| Item | Status |
|---|---|
| Interfaces to BoltzTraP and FermiSurfer, COHP, automatic k paths | Not yet written |
| Some solvers such as `Cluster-LNO` and EGAC | Partially only |
| Differences between upstream OpenMX and the GPU version | Not yet checked |
