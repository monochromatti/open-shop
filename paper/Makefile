TYPST ?= typst

.PHONY: all clean
all: paper.pdf

paper.pdf: paper.typ style.typ references.bib
	$(TYPST) compile --ignore-system-fonts paper.typ paper.pdf

clean:
	rm -f paper.pdf
