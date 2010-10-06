#include <stdio.h>
#include <stdlib.h>
#include <math.h>

enum bool { true, false };
int nextPrime(int, enum bool*) ;
void printPrimes(int, enum bool*) ;

int main(int argc, char** argv)
{
    int N, i;
    printf("Hi! Welcome to sequential sieve.\n");
    //printf("Enter N = ");
    if(argc != 2) {
        printf("Incorrect arguments idiot.\n");
        exit(1);
    }
    else
    {
        sscanf(argv[1], "%d", &N);
    }

    // --- begin

    enum bool *A = malloc(sizeof(enum bool)*N);

    A[0] = false;
    for(i = 1; i < N; i++)
    {
        A[i] = true;
    }

    int current = 2;

    while (current <= sqrt((double)N))
    {
        for(i = current; i <= N/current; i++)
        {
            A[i*current] = false;
        }
        current = nextPrime(current, A);
    }

    printPrimes(N, A);
}

int nextPrime(int c, enum bool*A) 
{
    int next = c+1; 
    while ( A[next] != true)
        next++;

    return next;
}

void printPrimes(int N, enum bool* A)
{
    int i;
    int nPrimes=0;
    for(i = 0; i < N; i++)
        if(A[i] == true)
        {
            printf("prime: %d\n", i);
            nPrimes++;
        }

    printf("So that's %d primes.\n", nPrimes); 
}
