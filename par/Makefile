CC= mpicc
#CFLAGS= -O3 
CFLAGS= -lbsponmpi
LFLAGS= -lm 
OBJIP= bspinprod.o bspedupack.o

all: ip 

ip: $(OBJIP)
	$(CC) $(CFLAGS) -o ip $(OBJIP) $(LFLAGS)

bspinprod.o: bspinprod.c bspedupack.h
	$(CC) $(CFLAGS) -c bspinprod.c

bspedupack.o: bspedupack.c bspedupack.h
	$(CC) $(CFLAGS) -c bspedupack.c

clean:
	rm -f *.o ip 
