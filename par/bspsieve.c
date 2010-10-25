#include "../bspedupack.h"
#include <limits.h>
#include <stdlib.h>

/* 
 * Author: Paul van der Walt
 * october 2010
 *
 * This program computes all prime numbers <= a given integer N.  
 */

// prototypes: 
ulong blockLow(int, int, ulong);
/*
ulong blockSize(ulong, ulong, ulong);
ulong blockHigh(ulong, ulong, ulong);
ulong blockOwner(ulong, ulong, ulong);
ulong globalIdx(ulong, ulong, ulong, ulong);
ulong localIdx(ulong, ulong, ulong, ulong);
*/

//globals:
int P; /* number of processors requested */ 
ulong N; /* requested max prime */

ulong findMinimum(int p, ulong* ks)
{

    ulong min = ULONG_MAX;
    ulong i;
    for(i = 0; i<p; i++)
    {
        min = MIN(ks[i], min); // range OK
    }

    return min;
}
ulong globalIdx(int p, int s, ulong n, ulong local)
{
    return local+blockLow(p,s,n);
}
ulong localIdx(int p, int s, ulong n, ulong global)
{
    return global-blockLow(p,s,n);
}
ulong blockLow(int p, int s, ulong n)
{
    // this means we've overflowed our max_value
    if(s*n < 0) 
         printf("Hm. s*n < 0: s=%lu, n=%lu, s*n=%lu\n", s,n,s*n);      
    return (s*n)/p; //implicit floor
}

ulong blockHigh(int p, int s, ulong n)
{
    return blockLow(p, s+1, n) - 1;
}

int blockOwner(int p, ulong index, ulong n)
{
    if(p*index < 0)
        printf("Hm. p*index < 0: p = %lu, i = %lu, p*i = %lu\n", p, index, p*index);
    return ((p*(index+1))-1)/n; //implicit floor
}

ulong blockSize(int p, int s, ulong n){
    /* Compute number of local components of processor s for vector
       of length n distributed over p processors with balanced
       block distribution (see paper). */

    return  blockLow(p,s+1,n)-blockLow(p,s,n) ; 

} /* end blockSize */

void bspmarkmultiples(int p, int s, ulong n, ulong k, ulong *x)
{
    // mark all multiples of k as non-prime in x
    /*
     * if (!(low_value % prime)) first = 0;
     * else first = prime - (low_value % prime);
     */

    ulong i;
    ulong first;

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

ulong nextPrime(int p, int s, ulong n, ulong k, ulong *x)
{
    // find minimal i s.t. i > k and i unmarked

    ulong newK = k+1;
    ulong local = MAX(
                    localIdx(p,s,n,newK),
                    0); // don't consider primes outside our range

    while(local < blockSize(p,s,n)-1 && x[local] == 0)
        local++;

    if(x[local] == 0)
    {
        return ULONG_MAX; // no primes for this processor. This is possible.
    }
    else
    {
        // if we get here we assume we found a prime!

        return globalIdx(p,s,n,local);
    }

} /* end nextPrime */

void bspsieve(){
    
    double alpha, time0, time1;
    ulong *x;  // local list of candidates
    ulong *ks; //place for proc0 to store intermediate k's
    ulong  
          n, 
          nl, 
          i, 
          iglob;
    int   s,
          p;
    ulong k;   // the current largest sure-prime

    n = N+1; // copy global N and increase by 1. (only proc 1 knows this)
             // this is so the maximum array idx == N
    
    bsp_begin(P);
    p= bsp_nprocs(); /* p = number of processors obtained */ 
    printf("Now we have %lu processors.\n", p);
    s= bsp_pid();    /* s = processor number */ 
    if (s==0){
        if(n<0)
            bsp_abort("Error in input: n is negative");
        ks = vecalloculi(p);
        if(ks == NULL)
            bsp_abort("Couldn't allocate ks[]!\n");
    }

    bsp_push_reg(&n,SZUL);
    bsp_sync();

    bsp_get(0,&n,0,&n,SZUL); //everyone reads N from proc 0
    bsp_sync();
    bsp_pop_reg(&n);

    printf("Sizeof(unsigned long int) == %zu\n", sizeof(unsigned int));
    nl= blockSize(p,s,n); // how big must s's block be?
    printf("P(%lu) tries to alloc vec of %lu ulongs = %lu Mb\n",s,  nl, nl*SZUL/1024/1024);
    x= vecalloculi(nl);
    if(x == NULL)
        bsp_abort("Couldn't allocate x[]!\n");

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

        bsp_push_reg(&k, SZUL);
        bsp_sync();

        if(s==0)
        {
            ks[0] = k; // my k
            for(i=1;i<p; i++)
            {
                bsp_get(i, &k, 0, &ks[i], SZUL);
            }
        }

        bsp_sync();

        if(s==0)
        {
            k = findMinimum(p,ks);
        }
        bsp_sync();

        //broadcast minimum 
        bsp_get(0,&k,0,&k,SZUL); 
        bsp_sync();

        bsp_pop_reg(&k);
    }

    // end work
    bsp_sync();  
    time1=bsp_time();

    ulong primes= 0;
    //printf("Processor %lu primes: \n", s); 
    for(i = 0; i < blockSize(p,s,n); i++)
        if( x[i] != 0)
            primes++;
     //       printf("  %lu is prime\n", globalIdx(p,s,n,i));
     printf("proc %lu finds %lu primes.\n", s, primes);

    fflush(stdout);
    if (s==0){
        printf("This took only %.6lf seconds.\n", time1-time0);
        fflush(stdout);
        vecfreeuli(ks);
    }

    vecfreeuli(x);
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
    sscanf(argv[1], "%lu", &N);

    printf("max prime requested = %lu\n", N);
    P = bsp_nprocs(); // maximum amount of procs

    if ( blockSize(P, 0, N) < sqrt(N))
        printf("WARNING: such a large P (%lu) with relatively small N (%lu) is inefficient. \n Choosing a lower P is recommended.\n\n", P, N);

    printf("Using %lu processors. \n", P);

    /* SPMD part */
    bspsieve();

    /* sequential part */
    exit(0);

} /* end main */
