# ab-initio-common-rust-substrate

**Can a small team, working with AI coding agents, build first-principles electronic-structure codes
that people can trust?** This repository holds the design and specifications for an experiment that
tries to answer that question, with Rust as the implementation language.

## The idea in one picture

```text
        LCAO (OpenMX-style)           PAW (plane waves)            FLAPW (all-electron)
   localized orbitals, H(R)/S(R)   accurate energies/forces    all-electron reference
   large systems, NEGF, ML-Ham.    e.g. ML-potential data      core levels, magnetism
             │                             │                            │
             └──────────── call shared libraries, keep own internals ───┘
                                           │
   radial solvers · spherical harmonics · XC · FFT · Ewald · symmetry & k-points · mixing · MD
                                           │
                        tenferro-rs: tensors and linear algebra on CPU/GPU
   ──────────────────────────────────────────────────────────────────────────────────────
   one verification system: specs per component · invariant tests · reference fixtures ·
   cross-method comparison · external all-electron data (ACWF: WIEN2k/FLEUR)
```

## Why

- **AI agents write code fast, but cannot vouch for it.** Correctness has to come from outside the agent:
  written contracts for every component and tests that decide whether each one is right.
- **Legacy codes are hard to verify piece by piece.** Agents can speed up an existing Fortran code
  (e.g. GPU patches to Quantum ESPRESSO), but a monolith can only be checked end to end against itself.
  A new code designed around component-level specifications can be checked part by part.
- **Permissive licensing matters.** Most established codes are GPL or proprietary. An MIT/Apache-2.0 code can be
  embedded in other software, ML pipelines and commercial tools.

## How correctness is established

1. **Specifications, not source code.** Each component (basis functions, integrals, Hamiltonian blocks, mixing,
   forces, ...) is specified with its inputs, outputs, conventions, thresholds and edge cases.
   Where the reference code is GPL, specifications are written under a clean-room procedure, and implementers
   never read the reference source.
2. **Four layers of tests.** Analytic solutions, physical invariants (Hermiticity, symmetries, sum rules,
   finite-difference forces), comparison with fixtures from reference codes, and SCF regression.
3. **Three methods check each other.** LCAO, PAW and FLAPW are compared on physical observables
   (equations of state, forces, bands), and against published all-electron reference data.

## What is shared, and what is not

The substrate is **a set of libraries, not a framework**. Methods share leaf components and outer contracts
(structure input, observables output, the verification harness), but each method owns its SCF loop, basis,
Hamiltonian storage and forces. Shared code is extracted only when real duplication appears.

## Status (November 2026)

| Part | Status |
|---|---|
| LCAO (OpenMX-style) | 18 specification chapters drafted and AI-reviewed; implementation not started |
| PAW, FLAPW | Scoping notes and prior-art survey (SIRIUS, FLEUR, dftatom, ...) |
| Rust implementation | Separate repositories, to be created |

The OpenMX-style method is used with the agreement of OpenMX's principal author; the implementation will carry
a different name and is not an official successor of OpenMX.

## Read more

- Book introduction: [`src/README.md`](src/README.md)
- Substrate overview: [`src/substrate/README.md`](src/substrate/README.md)
- Clean-room procedure: [`src/process/clean-room.md`](src/process/clean-room.md)
- Licensing of references and data: [`src/process/license.md`](src/process/license.md)
- LCAO specifications: [`src/methods/lcao-openmx/`](src/methods/lcao-openmx/README.md)

Build the book locally:

```sh
cargo install mdbook   # if not installed
mdbook serve --open
```

## Repository layout

- `src/process/` clean-room procedure, writing rules, licensing
- `src/conventions/` notation and units shared by all methods
- `src/testing/` test layers, fixtures, cross-method verification
- `src/substrate/` shared libraries and outer contracts; prior art
- `src/methods/` method families: `lcao-openmx/`, `paw/`, `flapw/`
- `tools/` scripts that check equations and conventions numerically
- `reviews/` clean-room review records

## License

- Text: [CC BY 4.0](LICENSE-CC-BY-4.0.txt)
- Code samples, scripts and test data: [MIT](LICENSE-MIT)

This repository contains no code from OpenMX or other GPL/proprietary codes.
