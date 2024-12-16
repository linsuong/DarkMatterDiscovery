
all: manual_6.1.pdf 

manual_6.1.pdf:manual_6.1.tex
	pdflatex  manual_6.1.tex
	bibtex manual_6.1
	pdflatex   manual_6.1.tex
	pdflatex   manual_6.1.tex
clean: 
	rm -f  *.dvi *.log *.aux *.toc *.bbl *.out *.blg .DS_Store