#include "../bspedupack.h"
#include <limits.h>

/* 
 * Author: Paul van der Walt
 * october 2010
 *
 * This program computes all prime numbers <= a given integer N.  
 */

// prototypes: 
int blockSize(int, int, int);
int blockLow(int, int, int);
int blockHigh(int, int, int);
int blockOwner(int, int, int);
int globalIdx(int, int, int, int);
int localIdx(int, int, int, int);
double bspip(int p, int s, int n, double *x, double *y);

//globals:
int P; /* number of processors requested */ 
int N; /* requested max prime */

int findMinimum(int p, int* ks)
{

    int min = INT_MAX;
    int i;
    for(i = 0; i<p; i++)
        min = MIN(ks[i], min);

    return min;
}
int globalIdx(int p, int s, int n, int local)
{
    return local+blockLow(p,s,n);
}
int localIdx(int p, int s, int n, int global)
{
    return global-blockLow(p,s,n);
}
int blockLow(int p, int s, int n)
{
    return (s*n)/p; //implicit floor
}

int blockHigh(int p, int s, int n)
{
    return blockLow(p, s+1, n) - 1;
}

int blockOwner(int p, int index, int n)
{
    return ((p*(index+1))-1)/n; //implicit floor
}

int blockSize(int p, int s, int n){
    /* Compute number of local components of processor s for vector
       of length n distributed over p processors with balanced
       block distribution (see paper). */

    return  blockLow(p,s+1,n)-blockLow(p,s,n) ; 

} /* end blockSize */

void bspmarkmultiples(int p, int s, int n, int k, int *x)
{
    // mark all multiples of k as non-prime in x
    /*
     * if (!(low_value % prime)) first = 0;
     * else first = prime - (low_value % prime);
     */

    int i;
    int first;

    if(k * k > blockLow(p,s,n)) // k is in our block, or beyond
        first = k * k - blockLow(p,s,n);
    else // k falls before our block; 
    {
        if(!(blockLow(p,s,n) % k)) // first element in our block is divisible by k
            first = 0;
        else
            first = k - (blockLow(p,s,n) % k); // start at first multiple of k in block
    }

    for (i=first; i < blockSize(p,s,n); i+= k)
    {
        x[i] = 0; //not a prime
    }


} /* end bspmarkmultiples */

int nextPrime(int p, int s, int n, int k, int *x)
{
    // find minimal i s.t. i > k and i unmarked

    int newK = k+1;
    int local = MAX(
                    localIdx(p,s,n,newK),
                    0); // don't consider primes outside our range

    while(local < blockSize(p,s,n)-1 && x[local] == 0)
        local++;

    if(x[local] == 0)
    {
        return INT_MAX; // no primes for this processor. This is possible.
    }
    else
    {
        // if we get here we assume we found a prime!

        return globalIdx(p,s,n,local);
    }

} /* end nextPrime */

void bspsieve(){
    
    double alpha, time0, time1;
    int *x;
    int *ks; //place for proc0 to store intermediate k's
    int p, s, n, nl, i, iglob;
    int k;   // the current largest sure-prime

    n = N+1; // copy global N and increase by 1. (only proc 1 knows this)
             // this is so the maximum array idx == N
    
    bsp_begin(P);
    p= bsp_nprocs(); /* p = number of processors obtained */ 
    s= bsp_pid();    /* s = processor number */ 
    if (s==0){
        if(n<0)
            bsp_abort("Error in input: n is negative");
        ks = vecalloci(p);
    }

    bsp_push_reg(&n,SZINT);
    bsp_sync();

    bsp_get(0,&n,0,&n,SZINT); //everyone reads N from proc 0
    bsp_sync();
    bsp_pop_reg(&n);

    nl= blockSize(p,s,n); // how big must s's block be?
    x= vecalloci(nl);
    for (i=0; i<nl; i++){
        // start by assuming everything is prime, except 1
        iglob= globalIdx(p,s,n,i);
        x[i]= iglob;
    }
    if(s==0)
        x[1]=0;
    bsp_sync(); 
    time0=bsp_time();
    k = 2;
    // begin work

    while( k*k <= n )
    {
        bspmarkmultiples(p,s,n,k,x);
        k = nextPrime(p,s,n,k,x);

        bsp_push_reg(&k, SZINT);
        bsp_sync();

        if(s==0)
        {
            ks[0] = k; // my k
            for(i=1;i<p; i++)
            {
                bsp_get(i, &k, 0, &ks[i], SZINT);
            }
        }

        bsp_sync();

        if(s==0)
        {
            k = findMinimum(p,ks);
        }
        bsp_sync();

        //broadcast minimum 
        bsp_get(0,&k,0,&k,SZINT); 
        bsp_sync();

        bsp_pop_reg(&k);
    }

    // end work
    bsp_sync();  
    time1=bsp_time();

    printf("Processor %d primes: \n", s); 
    for(i = 0; i < blockSize(p,s,n); i++)
        if( x[i] != 0)
            printf("  %d is prime\n", globalIdx(p,s,n,i));

    fflush(stdout);
    if (s==0){
        printf("This took only %.6lf seconds.\n", time1-time0);
        fflush(stdout);
        vecfreei(ks);
    }

    vecfreei(x);
    bsp_end();

} /* end bspsieve */

int main(int argc, char **argv){

    bsp_init(bspsieve, argc, argv);

    /* sequential part */
    if (argc != 2)
    {
        printf("Usage: %s N\n", argv[0]);
        bsp_abort("Incorrect invocation.\n");
    }
    sscanf(argv[1], "%d", &N);

    printf("max prime requested = %d\n", N);
    P = bsp_nprocs(); // maximum amount of procs

    if ( blockSize(P, 0, N) < sqrt(N))
        printf("WARNING: such a large P (%d) with relatively small N (%d) is inefficient. \n Choosing a lower P is recommended.\n\n", P, N);

    printf("Using %d processors. \n", P);

    /* SPMD part */
    bspsieve();

    /* sequential part */
    exit(0);

} /* end main */
