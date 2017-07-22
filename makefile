
run_flow:
	gcc -o flow -std=c99 flow.c
	./flow

run_rate: 
	gcc -o rate -std=c99 rate.c -lgsl -lgslcblas -lm
	./rate < matches.txt
