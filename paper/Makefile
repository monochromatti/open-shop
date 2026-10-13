LATEX ?= pdflatex

.PHONY: all clean
all: paper.pdf

paper.pdf: paper.tex
	mkdir -p .build
	$(LATEX) -interaction=nonstopmode -halt-on-error -output-directory=.build paper.tex
	$(LATEX) -interaction=nonstopmode -halt-on-error -output-directory=.build paper.tex
	cp .build/paper.pdf $@

clean:
	rm -rf .build
