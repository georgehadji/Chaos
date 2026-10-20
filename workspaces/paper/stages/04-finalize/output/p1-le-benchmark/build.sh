#!/bin/sh
# Rebuild main.tex and main.pdf from main.md (needs pandoc, pdflatex, bibtex). Hashes are wrapped with seqsplit so table cells break.
set -e
pandoc main.md -f markdown+pipe_tables -s --natbib --bibliography=refs.bib -H hdr.tex -V documentclass=article -V geometry:margin=2.5cm -V fontsize=11pt -o main.tex
python -c "import re;t=open('main.tex',encoding='utf-8').read();t=re.sub(r'(?<![0-9a-f])([0-9a-f]{40}|[0-9a-f]{64})(?![0-9a-f])',lambda m:'\\seqsplit{'+m.group(1)+'}',t);open('main.tex','w',encoding='utf-8').write(t)"
pdflatex -interaction=nonstopmode main.tex && bibtex main && pdflatex -interaction=nonstopmode main.tex && pdflatex -interaction=nonstopmode main.tex
