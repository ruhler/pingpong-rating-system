
run_flow:
	gcc -o flow -pedantic -std=c99 -ggdb -Wall -Werror flow.c
	./flow

run_rate: 
	gcc -o rate -std=c99 rate.c -lgsl -lgslcblas -lm
	./rate < matches.txt
