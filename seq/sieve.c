#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define ulong long

enum bool { true, false };
ulong nextPrime(ulong, enum bool*) ;
void printPrimes(ulong, enum bool*) ;

int main(int argc, char** argv)
{
    ulong N, i;
    printf("Hi! Welcome to sequential sieve.\n");
    if(argc != 2) {
        printf("Incorrect arguments idiot.\n");
        exit(1);
    }
    else
    {
        sscanf(argv[1], "%lu", &N);
        printf("you entered N== %lu\n", N);
    }

    // --- begin

    time_t t1,t2;
    time(&t1);

    printf("Trying to alloc %lu Mb... \n",N*sizeof(enum bool)/1024/1024); 
    enum bool *A = malloc(sizeof(enum bool)*N);
    if(A == NULL)
    { printf("Bad news, malloc failed. Try lower N.\n"); exit(1); }

    // assume all are prime to start with, 
    // except 0 and 1
    A[0] = false;
    A[1] = false;
    for(i = 2; i < N; i++)
    {
        A[i] = true;
    }

    // ... and start with k <- 2
    ulong current = 2;

    while (current*current <= N)
    {
        // start at current^2, optimisation
        for(i = current; i <= N/current; i++)
        {
            A[i*current] = false;
        }
        current = nextPrime(current, A);
    }

    printPrimes(N, A);
    time(&t2);
    printf("And it took: %lf sec\n", difftime(t2,t1)); 
    exit(0);
}

ulong nextPrime(ulong c, enum bool*A) 
{
    ulong next = c+1; 
    while ( A[next] != true)
        next++;

    return next;
}

void printPrimes(ulong N, enum bool* A)
{
    ulong i;
    ulong nPrimes=0;
    for(i = 0; i < N; i++)
        if(A[i] == true)
        {
            // do not print primes, just count them.
            nPrimes++;
        }

    printf("==> %lu primes.\n", nPrimes); 
}
