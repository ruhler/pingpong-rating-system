
.PHONY: all
all: pingpong-rating-system.pdf pingpong-rating-system

pingpong-rating-system.pdf: pingpong-rating-system.tex
	pdflatex pingpong-rating-system.tex
	pdflatex pingpong-rating-system.tex

pingpong-rating-system: pingpong-rating-system.c
	gcc -o $@ -std=c99 -Wall -Werror -Wshadow -ggdb $< -lm

