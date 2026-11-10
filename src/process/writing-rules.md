# Writing rules and source categories

## What may be written

- Definitions of physical quantities, equations, and their derivations (rewritten in our own notation).
- The idea of an algorithm (what is computed, in what order, and what it yields).
- Conventions as facts: units, signs, phases, ordering of indices and orbitals,
  file-format keywords and the meaning of columns.
- Numerical parameters: thresholds, grid sizes, cutoffs, with their values and meaning.
- Test inputs, expected values and tolerances.

## What must not be written (for GPL or proprietary reference codes)

- Fragments of the reference source code, or pseudocode paraphrasing it.
- Identifiers such as function, variable, or source-file names (file-format keywords are allowed).
- The way the implementation is assembled: function decomposition, loop structure, order of branches.
- Transcriptions of comments or documentation.
- Transcriptions of text, figures or tables from technical notes, papers or manuals.

## Public documents

- Rewrite equations, definitions and numbers in our own notation and cite the source
  (document and equation number).
- Quote the original only when necessary, briefly and with attribution.
- When a typo or inconsistency in a public document is found, record the verified result as a note
  in the [Public] section.

## Source-category tags

Each heading or paragraph carries one of the following tags.

| Tag | Meaning | Review |
|---|---|---|
| [Public] | Based on public documents | Check that sources are cited |
| [Convention] | Not in public documents; a fact about the behavior of an existing implementation | Reviewer confirms there is no code-derived expression |
| [Design] | A design decision made in this book | Not required |
| [Open] | Not yet investigated or confirmed | Not used as an implementation premise until resolved |

Qualifiers such as "(derived)" (a derivation by the spec writer) and "(checked)"
(numerically checked by a script in `tools/`) may follow a tag.

## No bitwise agreement

The goal of a reimplementation is agreement with the reference within the tolerances defined by the
tests, not bitwise agreement. Implementation details that only perturb numerical results slightly,
such as endpoint handling in interpolation, are in principle not written as [Convention].
Only what affects physical quantities or convergence (thresholds, cutoffs, regularization methods, ...)
is written.

## Writing math

- Inline math is written as `\\( ... \\)`; display math is enclosed by lines containing only `$$`.
- mdBook parses Markdown before passing it to MathJax, so `_`, `*`, `\,` etc. inside math need escaping.
  There is no need to do this by hand: `python3 tools/math_escape.py --fix src` converts automatically,
  and CI detects unconverted places with `--check`.
- Inside tables, do not use `|` in math; use `\lvert`, `\rvert`, `\lVert`, `\rVert`.

## Language

The book is written in English. Use straight ASCII punctuation in prose.
