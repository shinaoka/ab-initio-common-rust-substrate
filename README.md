# ab-initio-common-rust-substrate

Design and clean-room specifications for a common Rust substrate for ab initio electronic-structure methods:
an LCAO family in the style of [OpenMX](https://www.openmx-square.org/), plane-wave PAW, and FLAPW.
The book is written with [mdBook](https://rust-lang.github.io/mdBook/). Rust code lives in separate repositories.

## Build

```sh
cargo install mdbook   # if not installed
mdbook serve --open
```

## Layout

- `src/README.md` vision: why, what, and how correctness is established
- `src/process/` clean-room procedure, writing rules, licensing
- `src/conventions/` notation and units shared by all methods
- `src/testing/` test layers, fixtures, cross-method verification
- `src/substrate/` shared leaf libraries and outer contracts; prior art
- `src/methods/lcao-openmx/` LCAO (OpenMX-style) specifications: design, public documents, modules
- `src/methods/paw/`, `src/methods/flapw/` scoping notes
- `tools/` scripts that check equations and conventions numerically
- `reviews/` clean-room review records

## License

- Text: [CC BY 4.0](LICENSE-CC-BY-4.0.txt)
- Code samples, scripts and test data: [MIT](LICENSE-MIT)

This repository contains no code from OpenMX or other GPL/proprietary codes.
See `src/process/license.md` for the status of reference codes and data, including the OpenMX PAO/VPS databases.
