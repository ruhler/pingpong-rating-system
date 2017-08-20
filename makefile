
run_rate: 
	gcc -o rate -std=c99 rate.c
	./rate < matches.txt

run_fixed:
	gcc -o fixed -std=c99 -ggdb -Wall -Werror fixed.c
	./fixed

run_flow:
	gcc -o flow -std=c99 -ggdb -Wall -Werror flow.c
	./flow --test

run_prob: 
	gcc -o prob -std=c99 prob.c -lgsl -lgslcblas -lm
	./prob < matches.txt
