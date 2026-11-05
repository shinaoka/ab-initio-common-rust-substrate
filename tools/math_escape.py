# SPDX-License-Identifier: MIT
"""Escape TeX inside Markdown math so that mdBook's Markdown parser leaves it intact.

mdBook passes math to MathJax only after Markdown parsing.  Inside math regions
(\\( ... \\) inline, $$ ... $$ display) we therefore need
  * `_`  -> `\\_`   (otherwise pairs of underscores become <em>)
  * `*`  -> `\\*`
  * `\\,` `\\;` `\\!` `\\{` `\\}` `\\|` -> doubled backslash (Markdown eats a backslash before ASCII punctuation)

Usage:
  python3 tools/math_escape.py --check src   # exit 1 if anything needs escaping
  python3 tools/math_escape.py --fix src     # rewrite files in place
"""
import pathlib
import re
import sys

INLINE = re.compile(r"\\\\\((.+?)\\\\\)", re.S)
DISPLAY = re.compile(r"(^\$\$\n)(.*?)(^\$\$$)", re.S | re.M)
PUNCT = re.compile(r"(?<!\\)\\([,;!{}|])")
UNDERSCORE = re.compile(r"(?<!\\)_")
STAR = re.compile(r"(?<!\\)\*")


def escape(tex: str) -> str:
    tex = PUNCT.sub(r"\\\\\1", tex)
    tex = UNDERSCORE.sub(r"\\_", tex)
    tex = STAR.sub(r"\\*", tex)
    return tex


LEADING = re.compile(r"^(\s*)([-+=#>]|\d+\.)", re.M)


def escape_display(tex: str) -> str:
    # A display line starting with '-', '+', '=', '#', '>' or '1.' would be parsed
    # by Markdown as a list item / heading / quote; prefix an empty TeX group.
    return LEADING.sub(r"\1{}\2", escape(tex))


def process(text: str) -> str:
    text = DISPLAY.sub(lambda m: m.group(1) + escape_display(m.group(2)) + m.group(3), text)
    text = INLINE.sub(lambda m: "\\\\(" + escape(m.group(1)) + "\\\\)", text)
    return text


def main() -> int:
    mode, root = sys.argv[1], pathlib.Path(sys.argv[2])
    dirty = []
    for path in sorted(root.rglob("*.md")):
        old = path.read_text(encoding="utf-8")
        new = process(old)
        if new != old:
            dirty.append(path)
            if mode == "--fix":
                path.write_text(new, encoding="utf-8")
    for path in dirty:
        print(("fixed " if mode == "--fix" else "needs escaping: ") + str(path))
    return 1 if (dirty and mode == "--check") else 0


if __name__ == "__main__":
    sys.exit(main())
