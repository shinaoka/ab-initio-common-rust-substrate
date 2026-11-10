# LCAO family (OpenMX-style)

This part specifies an electronic-structure method in the style of
[OpenMX](https://www.openmx-square.org/): numerical pseudo-atomic orbitals (PAO) as the basis,
norm-conserving pseudopotentials (VPS), and real-space grids for the density.
The aim is a specification precise enough to **reimplement the method without reading the OpenMX source code**.

The implementation will not be called "OpenMX" and is not an official successor of OpenMX.

## Goals

- Rebuild the method as a set of components to which unit tests can be applied.
- Start from public technical notes, papers and the manual, and supply the conventions,
  numerical parameters and boundary conditions they leave out.
- Define, for each component, tests that decide correctness (analytic solutions, physical invariants,
  comparison with reference values).

## Role in the substrate

The family is centred on the real-space Hamiltonian and overlap H(R), S(R) in a
[sparse block format](design/sparse-block-format.md). This makes it:

- a data source for machine-learned Hamiltonians (e.g. DeepH), producing H(R)/S(R) for many small systems;
- a downstream solver for large-system H(R) predicted by such models (sparse diagonalization, O(N), NEGF, Green's functions);
- the home of OpenMX's established strengths: noncollinear magnetism, SOC, NEGF transport, O(N) methods, Wannier functions.

## Status

| Part | Status |
|---|---|
| Design, public-document index | First draft |
| 10 core module chapters (basis, two-center integrals, nonlocal pseudopotentials, neutral atom, grid/Poisson/XC, k-points, generalized eigenproblem, density matrix and mixing, total energy and forces, noncollinear) | First draft (AI-reviewed) |
| 6 advanced-feature chapters (O(N), DFT+U and magnetism, ESM/field/core level, EXX and dispersion, NEGF, geometry/MD) | First draft (AI-reviewed) |
| 2 analysis chapters | First draft (AI-reviewed) |

The main features of OpenMX 4.0 are covered. [Open] items in each chapter are to be settled against fixtures
([Gaps in public documents](sources/gaps.md)). Facts written as [Convention] were checked against
the GPU version of OpenMX by H. Kawai (commit 3871b1f) and may differ from upstream OpenMX.
Some equations and conventions are checked numerically by scripts in `tools/`.

## Reading order

1. [Module structure and data flow](design/architecture.md)
2. [H(R)/S(R) sparse block format](design/sparse-block-format.md)
3. Core modules, from [PAO basis](modules/pao-basis.md) to [noncollinear energy and forces](modules/noncollinear-energy-force.md)
4. Advanced features and analysis as needed
