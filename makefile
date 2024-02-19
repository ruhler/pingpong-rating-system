
pingpong-rating-system.pdf: pingpong-rating-system.tex
	pdflatex pingpong-rating-system.tex
	pdflatex pingpong-rating-system.tex

run_rate: 
	gcc -o rate -std=c99 -ggdb -Wall -Werror rate.c -lm
	./rate < matches.txt

