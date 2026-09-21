# Requires a TeX Live with pdflatex (Fedora: sudo dnf install texlive-scheme-medium).
# Overleaf also works: upload latex/ and ../FIGURES/ keeping the relative layout, or edit the
# \includegraphics paths.
paper.pdf: paper.tex
	pdflatex -interaction=nonstopmode paper.tex
	pdflatex -interaction=nonstopmode paper.tex   # 2nd pass: xltabular widths, references

paper.tex: ../PAPER_MANUSCRIPT_V2.md md2tex.py preamble.tex
	python3 md2tex.py
