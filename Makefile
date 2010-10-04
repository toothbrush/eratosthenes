#CFLAGS=-lm

all: sieve

sieve: sieve.o
	gcc -lm sieve.o -o sieve

sieve.o: sieve.c
	gcc -c sieve.c

clean:
	rm -v sieve{,.o}
