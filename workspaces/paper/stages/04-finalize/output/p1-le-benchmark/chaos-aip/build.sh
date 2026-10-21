#!/bin/sh
# Rebuild main.pdf from main.tex (REVTeX 4.2, AIP Chaos style; needs pdflatex and bibtex).
set -e
pdflatex -interaction=nonstopmode main.tex && bibtex main && pdflatex -interaction=nonstopmode main.tex && pdflatex -interaction=nonstopmode main.tex
# build leftovers; mainNotes.bib breaks the repo naming rule
rm -f main.aux main.blg main.log main.out mainNotes.bib
