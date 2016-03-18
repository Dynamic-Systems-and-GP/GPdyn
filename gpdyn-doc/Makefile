all:	readme pdf

readme:
	pandoc -s --mathjax -S -i manual.tex  -o ../README.rst

pdf:
	pdflatex manual.tex
	rm manual.log manual.aux
