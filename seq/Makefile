#CFLAGS=-lm

all: sieve pretty

pretty:
	lgrind -i sieve.c > sieve.lg
sieve: sieve.o
	gcc -lm sieve.o -o sieve

sieve.o: sieve.c
	gcc -c sieve.c

clean:
	rm -f sieve{,.o}
