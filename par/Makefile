CC= mpicc
CFLAGS= -g -Wall
LFLAGS= -lm -lbsponmpi #-Wl,-rpath -Wl,LIBDIR
OBJSIEVE= bspsieve.o bspedupack.o

all: sieve pretty

pretty: bspsieve.c
	lgrind -i bspsieve.c > bspsieve.lg

sieve: $(OBJSIEVE)
	$(CC) $(CFLAGS) -o sieve $(OBJSIEVE) $(LFLAGS)

bspsieve.o: bspsieve.c ../bspedupack.h
	$(CC) $(CFLAGS) -c bspsieve.c

bspedupack.o: ../bspedupack.c ../bspedupack.h
	$(CC) $(CFLAGS) -c ../bspedupack.c

clean:
	rm -f *.o sieve
