
run_fixed:
	gcc -o fixed -std=c99 -ggdb -Wall -Werror fixed.c
	./fixed

run_flow:
	gcc -o flow -std=c99 -ggdb -Wall -Werror flow.c
	./flow --test

run_rate: 
	gcc -o rate -std=c99 rate.c -lgsl -lgslcblas -lm
	./rate < matches.txt
