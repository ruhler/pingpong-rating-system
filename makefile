
run_rate: 
	gcc -o rate -std=c99 rate.c -lgsl -lgslcblas -lm
	./rate
