# Setup: LaTeX (paper/04 only)

LaTeX typesets the final manuscript. Optional: without it, paper/04 still produces the LaTeX sources.

## Install

Windows: MiKTeX (https://miktex.org/) or TeX Live (https://tug.org/texlive/). Allow on-the-fly package installs in MiKTeX.

## Verify

`pdflatex --version` and `bibtex --version` print versions.

## How this repo uses it

paper/04 converts the revised draft to the journal template and builds the PDF in `workspaces/paper/stages/04-finalize/output/[paper]/`.
