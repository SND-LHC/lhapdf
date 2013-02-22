LATEXCMD := pdflatex -output-directory=output lhapdf6-manual.tex

lhapdf6-manual.pdf: $(shell ls *.tex)
	$(LATEXCMD) && $(LATEXCMD) && cp output/lhapdf6-manual.pdf .

clean:
	rm -rf output
