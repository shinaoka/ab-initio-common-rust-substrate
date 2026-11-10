# Module structure and data flow [Design]

## Policy

- No module holds global state; inputs and outputs are made explicit through types.
- The central object passed between modules is the [H(R)/S(R) sparse block format](sparse-block-format.md).
- Each module is made testable in isolation ([test hierarchy](../../../testing/strategy.md)).

## Data flow

```text
Input (structure, basis specification, calculation conditions)
  │
  ├─► Basis: reading PAO and VPS, spherical Bessel transform of radial functions
  │
  ├─► Two-center integrals: S(R), T(R), nonlocal term, neutral-atom term
  │        │
  │        ▼
  │   Sparse block format H₀(R), S(R) (the part determined by structure and basis only)
  │
  └─► SCF loop
         ρ(r) ──► Grid processing (Poisson, XC, matrix elements) ──► H[ρ](R)
                                                              │
         ┌────────────────────────────────────────────────────┘
         ▼
      H(k), S(k) per k point ──► generalized eigenvalue ──► density matrix ρ(R) ──► ρ(r)
         │                                                        │
         └──────────────── charge mixing ◄────────────────────────┘
  │
  └─► Total energy, forces, stress, analysis (post-processing taking H(R), S(R) as input)
```

## Reuse as components

| Component | Input | Output | Anticipated external use |
|---|---|---|---|
| Basis and two-center integrals | Structure, basis specification | S(R), T(R), nonlocal term | Machine-learning of Hamiltonians, tight-binding models |
| Evaluation of H[ρ] | ρ, structure | H(R) | Automatic differentiation, coupling with machine learning |
| Generalized eigenvalue and O(N) methods | H(k), S(k) or H(R), S(R) | Eigenvalues and eigenvectors, density matrix | Other localized-basis codes |
| Post-processing | H(R), S(R) | Bands, exchange interactions, transport, etc. | Other localized-basis codes |
