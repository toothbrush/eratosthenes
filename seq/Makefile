all: sieve pretty

pretty:
	lgrind -i sieve.c > sieve.lg

sieve: sieve.o
	gcc -lm sieve.o -o sieve -Wall

sieve.o: sieve.c
	gcc -c sieve.c -Wall

clean:
	rm -f sieve{,.o}
