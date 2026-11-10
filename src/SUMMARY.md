# Summary

[Introduction](README.md)

# Process

- [Clean-room procedure](process/clean-room.md)
- [Writing rules and source categories](process/writing-rules.md)
- [Licensing](process/license.md)

# Common

- [Notation, units and terminology](conventions/notation.md)
- [Testing strategy](testing/strategy.md)
- [Substrate overview](substrate/README.md)
- [Prior art](substrate/prior-art.md)

# LCAO family (OpenMX-style)

- [Overview](methods/lcao-openmx/README.md)
- [Public documents](methods/lcao-openmx/sources/index.md)
- [Gaps in public documents](methods/lcao-openmx/sources/gaps.md)
- [Module structure and data flow](methods/lcao-openmx/design/architecture.md)
- [H(R)/S(R) sparse block format](methods/lcao-openmx/design/sparse-block-format.md)
- [Parallelization](methods/lcao-openmx/design/parallelization.md)
- [GPU strategy and porting policy](methods/lcao-openmx/design/gpu-strategy.md)
- [Core modules]()
  - [PAO basis and real spherical harmonics](methods/lcao-openmx/modules/pao-basis.md)
  - [Two-center integrals: overlap and kinetic energy](methods/lcao-openmx/modules/two-center-overlap-kinetic.md)
  - [Nonlocal pseudopotentials](methods/lcao-openmx/modules/nonlocal-pseudopotential.md)
  - [Neutral-atom potential and projector expansion](methods/lcao-openmx/modules/neutral-atom.md)
  - [Real-space grid, Poisson and exchange-correlation](methods/lcao-openmx/modules/grid-poisson-xc.md)
  - [k-points and Bloch sums](methods/lcao-openmx/modules/kpoints-bloch.md)
  - [Generalized eigenvalue problem and regularization of S](methods/lcao-openmx/modules/generalized-eigen.md)
  - [Density matrix and charge mixing](methods/lcao-openmx/modules/density-mixing.md)
  - [Total energy and forces](methods/lcao-openmx/modules/energy-force.md)
  - [Noncollinear and spin-orbit energy, forces and stress](methods/lcao-openmx/modules/noncollinear-energy-force.md)
- [Advanced features]()
  - [O(N) methods](methods/lcao-openmx/modules/on-methods.md)
  - [DFT+U, constrained DFT and magnetic features](methods/lcao-openmx/modules/dftu-magnetism.md)
  - [ESM, electric field, charge doping, virtual atoms, core levels](methods/lcao-openmx/modules/esm-field-corehole.md)
  - [Exact exchange, hybrid functionals and dispersion corrections](methods/lcao-openmx/modules/exx-vdw.md)
  - [Electron transport by NEGF](methods/lcao-openmx/modules/negf.md)
  - [Geometry optimization, molecular dynamics and NEB](methods/lcao-openmx/modules/geometry-md.md)
  - [Mixed-precision GPU numerics](methods/lcao-openmx/modules/gpu-mixed-precision.md)
- [Analysis]()
  - [Bands, DOS, charge analysis, exchange couplings, Berry phase, optics](methods/lcao-openmx/modules/analysis.md)
  - [Wannier functions, band unfolding, localized-orbital analysis](methods/lcao-openmx/modules/wannier-unfolding.md)

# PAW family

- [Overview](methods/paw/README.md)

# FLAPW family

- [Overview](methods/flapw/README.md)
