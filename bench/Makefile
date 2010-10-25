CC= mpicc
#CFLAGS= -O3 
CFLAGS= 
LFLAGS= -lbsponmpi -lm #-Wl,-rpath -Wl,LIBDIR
#OBJIP= bspinprod.o bspedupack.o
OBJBENCH= bspbench.o bspedupack.o

all: bench

bench: $(OBJBENCH)
	$(CC) $(CFLAGS) -o bench $(OBJBENCH) $(LFLAGS)

bspbench.o: bspbench.c ../bspedupack.h
	$(CC) $(CFLAGS) -c bspbench.c

bspedupack.o: ../bspedupack.c ../bspedupack.h
	$(CC) $(CFLAGS) -c ../bspedupack.c

clean:
	rm -f *.o bench
