
rate.pdf: rate.tex
	pdflatex rate.tex
	pdflatex rate.tex

run_rate: 
	gcc -o rate -std=c99 -ggdb -Wall -Werror rate.c -lm
	./rate < matches.txt

