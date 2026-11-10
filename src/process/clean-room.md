# Clean-room procedure

This procedure applies whenever a specification draws on a **reference code** whose license
does not allow its code to be reused under MIT/Apache-2.0 (GPL codes such as OpenMX,
Quantum ESPRESSO, ABINIT, GPAW, Elk; proprietary codes such as WIEN2k).
For permissively licensed references (e.g. SIRIUS, FLEUR, dftatom; see
[Prior art](../substrate/prior-art.md)) the procedure is not required, but their license
and copyright notices must be respected and recorded in the chapter that uses them.

## Roles

| Role | May read | Must not read | Output |
|---|---|---|---|
| Spec writer | The reference code, public documents, execution results | — | [Convention] sections of this book, fixtures |
| Reviewer | This book, the reference code | — | Confirmation that no code-derived expression has leaked into the book |
| Implementer | This book, public documents (papers, technical notes, manuals), fixtures | The reference code (including derived versions) | The new implementation |

- The spec writer and the implementer are different people (or different agent sessions).
  If the same person does both, it is hard to show that the implementation does not depend on the code.
- When agents are used, the implementer's session has no checkout of the reference code, and its
  network access is limited to the sites of public documents. Session logs are kept.
- The reviewer may be an AI agent. In that case it runs in a fresh session, separate from the one
  that wrote the chapter, compares the chapter with the reference code and public documents,
  and records the result.
- If the implementer finds the specification insufficient, they open an issue in this repository and
  the spec writer updates the book. The implementer does not inspect the reference code.

## Flow

1. The spec writer organizes equations and methods from public documents ([Public]).
2. The spec writer determines conventions and numerical parameters missing from the public
   documents by examining the reference implementation, and records them as facts ([Convention]).
3. The spec writer runs the reference implementation to produce fixtures (inputs and expected values).
4. The reviewer confirms that the [Convention] text contains no code-derived expression
   (code fragments, identifiers, transcribed control flow).
5. The implementer implements from the book, public documents and fixtures only, and makes the tests pass.

## Records

- Each version of the book is tracked in this repository's git history.
- Reviews are recorded in `reviews/` and in pull-request reviews.
- Implementation repositories record which version of the book they were built from.

## Notes

- The legal standing of AI-assisted clean-room implementation is not settled. This procedure leaves
  evidence against dependence on the reference code; it is not a legal guarantee.
- Language models have been trained on public scientific codes, so a perfect clean room is impossible.
  The aim is that implementations are designed from the specification and do not copy or adapt
  the structure of the reference code.
- Copyright holders of reference codes are consulted on the approach where practical.
  For OpenMX, see [Licensing](license.md).
