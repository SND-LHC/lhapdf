LATEXCMD := pdflatex -output-directory=output lhapdf6-manual.tex

lhapdf6-manual.pdf: $(shell ls *.tex)
	test -e output || mkdir output
	$(LATEXCMD) && $(LATEXCMD) && cp output/lhapdf6-manual.pdf .

clean:
	@rm -rf output auto lhapdf6-manual.pdf
