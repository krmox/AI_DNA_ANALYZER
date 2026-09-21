# LaTeX build of PAPER_MANUSCRIPT_V2.md

`paper.tex` is generated: `python3 md2tex.py` reads `../PAPER_MANUSCRIPT_V2.md` (the single source of truth,
never edit `paper.tex` by hand) and `preamble.tex` (layout). Build: `make` (pdflatex, two passes).

Status: **not compiled in the environment that produced it (no TeX installed).** Structural checks passed
(brace and environment balance, table column counts, no unmapped Unicode, figure files exist). Expect to fix
small warnings on the first real compile (overfull boxes in wide tables are the likely ones).

Yellow boxes are unresolved placeholders (`AUTHOR TO COMPLETE`, `CITATION REQUIRED`, `REFERENCE NEEDED`).
This is a neutral single-column article layout, not a journal template; a journal submission needs that
journal's own class, and md2tex.py's output would have to be adapted to it.
