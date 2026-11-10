# Substrate overview [Design]

## Principle: libraries, not a framework

All method families live in the same project, use the same tools and meet the same verification
standards, but **each keeps its own internal implementation and interfaces**.
The substrate provides leaf libraries that methods *call*, and outer contracts that methods *satisfy*.
It does not provide a common SCF driver or a common Hamiltonian abstraction.

Shared components are extracted only when real duplication appears between two method families,
not in advance.

```text
 ┌───────────────┐   ┌───────────────┐   ┌───────────────┐
 │ LCAO (OpenMX) │   │      PAW      │   │     FLAPW     │   each owns: SCF loop, basis,
 │               │   │               │   │               │   Hamiltonian storage, eigensolver
 │               │   │               │   │               │   usage, forces and stress
 └──────┬────────┘   └──────┬────────┘   └──────┬────────┘
        │  calls             │  calls             │  calls
 ┌──────┴────────────────────┴────────────────────┴───────┐
 │ leaf libraries: radial solvers, Ylm/Gaunt, XC, FFT,     │
 │ Ewald, symmetry and k-points, mixers, optimizers and MD │
 └─────────────────────────────┬───────────────────────────┘
                        tenferro-rs (tensors, linear algebra, CPU/GPU)
 ─────────────────────────────────────────────────────────────
 outer contracts: structure input · observables output · verification harness
```

## Candidate leaf libraries

| Component | Content | Used by |
|---|---|---|
| Structure and symmetry | Lattice, atoms, units, space group, irreducible k-points, band paths | all |
| Uniform real-space grid and FFT | Grid definition, CPU/GPU FFT, interpolation | all (LCAO: Hartree; PAW: everything; FLAPW: interstitial) |
| Hartree and XC on grids | Poisson solvers, LDA/GGA/meta-GGA, spin and noncollinear handling | all |
| Ewald sum | Ion-ion energy, forces and stress | all |
| Radial atomic data | Radial grids, spherical Bessel transforms, real spherical harmonics, Gaunt coefficients, Lebedev angular grids | all |
| Radial atomic solvers | Schrödinger, scalar-relativistic and Dirac radial equations; atomic SCF | FLAPW (core, radial functions), PAW (dataset generation), LCAO (PAO/VPS generation) |
| Separable nonlocal operators | Projectors \\( \lvert p \rangle D \langle p \rvert \\) and their coefficients | LCAO, PAW |
| On-site quantities | Occupation matrices, DFT+U, on-site SOC | all |
| SCF numerics | Pulay/Broyden/Kerker mixing, smearing, Fermi-level search, convergence checks | all |
| Eigensolver building blocks | Generalized eigenproblems, Rayleigh-Ritz, Cholesky orthonormalization, Davidson/LOBPCG | all, used differently |
| Moving atoms | BFGS/FIRE, MD integrators, thermostats, cell relaxation | all |

## Outer contracts

- **Structure input:** a single description of lattice, atoms, species and calculation conditions.
- **Observables output:** total energy and its components, forces, stress, eigenvalues and band structures,
  densities of states, magnetic moments, in internal units with documented conventions.
- **Verification harness:** reference-data management, tolerance specifications, invariant tests,
  cross-method comparison of observables, CI.
- **Development conventions:** repository layout (Cargo workspaces), agent instructions, test style.

## What is not shared

| | LCAO (OpenMX-style) | PAW | FLAPW |
|---|---|---|---|
| Wavefunctions | Atomic-orbital coefficients | Plane-wave coefficients | LAPW + local-orbital coefficients |
| Hamiltonian | Explicit sparse H(R), S(R) | Matrix-free \\( H\lvert\psi\rangle \\) via FFT | Explicit dense H, S |
| Characteristic terms | Neutral-atom potential, O(N) methods, NEGF | Compensation charges, on-site corrections | Matching conditions, local orbitals, linearization energies |
| Basis corrections | Pulay forces | Pulay stress | Pulay and core corrections to forces |

## Interfaces between methods

Methods are compared through observables only (see [Testing strategy](../testing/strategy.md)).
Within a method family, internal contracts are designed freely; for the LCAO family the central one is the
[H(R)/S(R) sparse block format](../methods/lcao-openmx/design/sparse-block-format.md).
