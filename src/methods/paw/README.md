# PAW family (plane waves) [Open]

Scoping only; no specification chapters yet.

## Purpose

A GPU-native plane-wave PAW implementation under MIT/Apache-2.0. The main demand is high-throughput,
accurate total energies, forces and stress, e.g. training data for machine-learned interatomic potentials,
a role currently dominated by proprietary codes.

## Why it fits the verification experiment

- External all-electron references exist that are independent of any single code
  (ACWF verification data, WIEN2k/FLEUR), so tests can check physical correctness, not only agreement with one code.
- With shared PAW datasets, components (projectors, compensation charges, on-site terms) can be compared one by one
  with other codes.

## Components specific to this family

- Plane-wave basis, cutoff sphere, matrix-free \\( H\lvert\psi\rangle \\) and \\( S\lvert\psi\rangle \\) via FFT.
- PAW datasets: reading (PAW-XML, UPF), partial waves, projectors, compensation charges.
- On-site densities and energies; on-site XC on Lebedev grids.
- Iterative eigensolvers (Davidson, LOBPCG); Pulay stress.
- Later: ultrasoft pseudopotentials, hybrid functionals with PAW.

## References to study

- SIRIUS (BSD-3-Clause): primary algorithm and design reference; also a component-level oracle.
- NWChem `nwpw/paw` (ECL-2.0): readable PAW implementation.
- DFTK.jl (MIT): plane-wave SCF and differentiability design.
- GPL codes (QE, ABINIT, GPAW, libPAW): clean-room references only.

See [Prior art](../../substrate/prior-art.md) for licenses.

## [Open]

- Source of PAW datasets: PseudoDojo/JTH redistribution terms are unknown; generating our own datasets
  with a radial solver (e.g. a port of dftatom) is an alternative.
- First milestone (e.g. equations of state for a subset of the ACWF set).
