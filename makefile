
run_rate: 
	gcc -o rate -std=c99 -ggdb -Wall -Werror -lgsl -lgslcblas -lm rate.c flow.c prob.c fixed.c
	./rate < matches.txt

