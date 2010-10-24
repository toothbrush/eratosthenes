CC= mpicc
#CFLAGS= -O3 
CFLAGS= -lbsponmpi
LFLAGS= -lm #-Wl,-rpath -Wl,LIBDIR
#OBJIP= bspinprod.o bspedupack.o
OBJSIEVE= bspsieve.o bspedupack.o

all: sieve

#ip: $(OBJIP)
#	$(CC) $(CFLAGS) -o ip $(OBJIP) $(LFLAGS)

sieve: $(OBJSIEVE)
	$(CC) $(CFLAGS) -o sieve $(OBJSIEVE) $(LFLAGS)

bspsieve.o: bspsieve.c ../bspedupack.h
	$(CC) $(CFLAGS) -c bspsieve.c

#bspinprod.o: bspinprod.c bspedupack.h
#	$(CC) $(CFLAGS) -c bspinprod.c

bspedupack.o: ../bspedupack.c ../bspedupack.h
	$(CC) $(CFLAGS) -c ../bspedupack.c

clean:
	rm -f *.o ip sieve
