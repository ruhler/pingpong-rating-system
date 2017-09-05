
rate.pdf: rate.tex
	pdflatex rate.tex
	pdflatex rate.tex

run_rate: 
	gcc -o rate -std=c99 -ggdb -Wall -Werror -lgsl -lgslcblas -lm rate.c prob2.c
	./rate < matches.txt

