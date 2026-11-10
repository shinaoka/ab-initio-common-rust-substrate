# Introduction

This book describes a **common Rust substrate for ab initio electronic-structure methods**
and the implementation specifications of the methods built on it. It serves two readers at once:

- **People** who want to understand the project, judge whether it is sound, and join it.
- **Implementers, human or AI agent,** who need specifications precise enough to write
  and test code without reading the source of existing GPL or proprietary codes.

The Rust code lives in separate repositories. This repository holds only the design,
the specifications, and the verification plan.

## Why

AI coding agents can now write large amounts of scientific code quickly. What they cannot
do on their own is guarantee that the code is *correct*. The question this project asks is:

> Can a small group, working with AI agents, build and maintain trustworthy
> first-principles codes, if correctness is pinned down by specifications and tests
> at the level of individual components?

If the answer is yes, implementing a new method or maintaining an existing one becomes far
cheaper, and small research groups can own reliable electronic-structure codes. We will publish
what we learn, including what fails.

## What

Three method families share one substrate, but each keeps its own internals and interfaces:

| Method family | Basis | Role | Status |
|---|---|---|---|
| [LCAO (OpenMX-style)](methods/lcao-openmx/README.md) | Numerical pseudo-atomic orbitals, norm-conserving pseudopotentials | H(R)/S(R)-centric electronic structure; data source and downstream solver for ML Hamiltonians; large systems, O(N), NEGF | Specification drafted (18 chapters) |
| [PAW](methods/paw/README.md) | Plane waves + projector augmentation | High-throughput, accurate total energies and forces (e.g. training data for ML potentials) | Scoping |
| [FLAPW](methods/flapw/README.md) | Linearized augmented plane waves, all-electron | All-electron reference within the same framework; core-level, hyperfine, magnetism | Scoping |

The substrate is **a set of libraries, not a framework**: shared leaf components
(tensor and linear algebra via tenferro-rs, FFT, radial solvers, spherical harmonics, XC,
symmetry, k-points, mixing) plus shared outer contracts (structure input, physical observables,
verification harness). Each method owns its SCF loop, basis representation, Hamiltonian
storage, and forces. See [Substrate overview](substrate/README.md).

## How correctness is established

1. **Component-level specifications.** Each component has a written contract: inputs,
   outputs, conventions, thresholds, edge cases.
2. **Four layers of tests** ([Testing strategy](testing/strategy.md)): analytic solutions,
   physical invariants (Hermiticity, translational and rotational invariance, sum rules,
   finite-difference forces), comparison with fixtures from reference codes, and SCF regression.
3. **Cross-method verification.** The three families are compared at the level of physical
   observables (equations of state, forces, band structures), and against external all-electron
   references such as the ACWF verification data set (WIEN2k/FLEUR).
4. **Clean-room discipline** where reference codes are GPL or proprietary
   ([Clean-room procedure](process/clean-room.md)).

## Licensing

- Text: CC BY 4.0. Code samples, scripts, test data: MIT.
- Implementations built from these specifications are intended to be MIT/Apache-2.0.
- See [Licensing](process/license.md) for the status of each reference code and data set,
  including OpenMX's PAO/VPS databases.

## How to read

Statements in method chapters carry a source-category tag
([Writing rules](process/writing-rules.md)):

- **[Public]** Based on public documents (papers, technical notes, manuals).
- **[Convention]** Not in public documents; a fact about the behavior of an existing implementation
  (units, ordering, thresholds, ...). Subject to clean-room review.
- **[Design]** A design decision made in this book, independent of any existing implementation.
- **[Open]** Not yet investigated or confirmed.

## Status and history

The LCAO chapters were first written in Japanese as "OpenMX implementation notes" and reviewed by
AI reviewers in that form (records in `reviews/`). They were translated into English when the
repository was reorganized into this substrate (November 2026).
