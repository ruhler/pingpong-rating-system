
run_rate: 
	gcc -o rate -std=c99 -ggdb -Wall -Werror rate.c flow.c
	./rate < matches.txt

run_fixed:
	gcc -o fixed -std=c99 -ggdb -Wall -Werror fixed.c
	./fixed

run_prob: 
	gcc -o prob -std=c99 prob.c -lgsl -lgslcblas -lm
	./prob < matches.txt
