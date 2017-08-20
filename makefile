
run_rate: 
	gcc -o rate -std=c99 -ggdb -Wall -Werror -lgsl -lgslcblas -lm rate.c flow.c prob.c
	./rate < matches.txt

run_fixed:
	gcc -o fixed -std=c99 -ggdb -Wall -Werror fixed.c
	./fixed

