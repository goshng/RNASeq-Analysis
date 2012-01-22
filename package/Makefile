LATEX=../../latex
export BIBINPUTS=.:${LATEX}/bib
export BSTINPUTS=.:${LATEX}/bst:
export TEXINPUTS=.:${LATEX}/sty:

OUTEXTS=aux bbl blg log pdf dvi lof lot toc

.SUFFIXES: .tex .pdf .Rnw

DOCS = supplementary manuscript 
#DOCS = supplementary

all: ${DOCS:%=%.pdf}

.Rnw.tex: 
	R CMD Sweave $*.Rnw

.tex.pdf: 
	R CMD pdflatex $*
	bibtex $*
	cp $*.tex a.tex
	R CMD pdflatex $*
	R CMD pdflatex $*
	open $*.pdf

clean:
	rm -f $(foreach d,${DOCS},${OUTEXTS:%=$d.%})
	rm -f ${DOCS:%=%*.pdf}
