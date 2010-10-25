#include "../bspedupack.h"
#include <limits.h>

/* 
 * Author: Paul van der Walt
 * october 2010
 *
 * This program computes all prime numbers <= a given integer N.  
 */

// prototypes: 
ulong blockSize(ulong, ulong, ulong);
ulong blockLow(ulong, ulong, ulong);
ulong blockHigh(ulong, ulong, ulong);
ulong blockOwner(ulong, ulong, ulong);
ulong globalIdx(ulong, ulong, ulong, ulong);
ulong localIdx(ulong, ulong, ulong, ulong);

//globals:
ulong P; /* number of processors requested */ 
ulong N; /* requested max prime */

ulong findMinimum(ulong p, ulong* ks)
{

    ulong min = ULONG_MAX;
    ulong i;
    for(i = 0; i<p; i++)
        min = MIN(ks[i], min);

    return min;
}
ulong globalIdx(ulong p, ulong s, ulong n, ulong local)
{
    return local+blockLow(p,s,n);
}
ulong localIdx(ulong p, ulong s, ulong n, ulong global)
{
    return global-blockLow(p,s,n);
}
ulong blockLow(ulong p, ulong s, ulong n)
{
    // this means we've overflowed our max_value
    if(s*n < 0) 
         printf("Hm. s*n < 0: s=%ld, n=%ld, s*n=%ld\n", s,n,s*n);      
    return (s*n)/p; //implicit floor
}

ulong blockHigh(ulong p, ulong s, ulong n)
{
    return blockLow(p, s+1, n) - 1;
}

ulong blockOwner(ulong p, ulong index, ulong n)
{
    return ((p*(index+1))-1)/n; //implicit floor
}

ulong blockSize(ulong p, ulong s, ulong n){
    /* Compute number of local components of processor s for vector
       of length n distributed over p processors with balanced
       block distribution (see paper). */

    return  blockLow(p,s+1,n)-blockLow(p,s,n) ; 

} /* end blockSize */

void bspmarkmultiples(ulong p, ulong s, ulong n, ulong k, ulong *x)
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

ulong nextPrime(ulong p, ulong s, ulong n, ulong k, ulong *x)
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
    ulong *x;
    ulong *ks; //place for proc0 to store intermediate k's
    ulong p, s, n, nl, i, iglob;
    ulong k;   // the current largest sure-prime

    n = N+1; // copy global N and increase by 1. (only proc 1 knows this)
             // this is so the maximum array idx == N
    
    bsp_begin(P);
    p= bsp_nprocs(); /* p = number of processors obtained */ 
    s= bsp_pid();    /* s = processor number */ 
    if (s==0){
        if(n<0)
            bsp_abort("Error in input: n is negative");
        ks = vecalloculi(p);
    }

    bsp_push_reg(&n,sizeof(ulong));
    bsp_sync();

    bsp_get(0,&n,0,&n,sizeof(ulong)); //everyone reads N from proc 0
    bsp_sync();
    bsp_pop_reg(&n);

    nl= blockSize(p,s,n); // how big must s's block be?
    printf("P(%ld) tries to alloc vec of %ld ulongs = %ld bytes\n",s,  nl, nl*sizeof(ulong));
    x= vecalloculi(nl);

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

        bsp_push_reg(&k, sizeof(ulong));
        bsp_sync();

        if(s==0)
        {
            ks[0] = k; // my k
            for(i=1;i<p; i++)
            {
                bsp_get(i, &k, 0, &ks[i], sizeof(ulong));
            }
        }

        bsp_sync();

        if(s==0)
        {
            k = findMinimum(p,ks);
        }
        bsp_sync();

        //broadcast minimum 
        bsp_get(0,&k,0,&k,sizeof(ulong)); 
        bsp_sync();

        bsp_pop_reg(&k);
    }

    // end work
    bsp_sync();  
    time1=bsp_time();

    ulong primes= 0;
    //printf("Processor %ld primes: \n", s); 
    for(i = 0; i < blockSize(p,s,n); i++)
        if( x[i] != 0)
            primes++;
     //       printf("  %ld is prime\n", globalIdx(p,s,n,i));
     printf("proc %ld finds %ld primes.\n", s, primes);

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
    sscanf(argv[1], "%ld", &N);

    printf("max prime requested = %ld\n", N);
    P = bsp_nprocs(); // maximum amount of procs

    if ( blockSize(P, 0, N) < sqrt(N))
        printf("WARNING: such a large P (%ld) with relatively small N (%ld) is inefficient. \n Choosing a lower P is recommended.\n\n", P, N);

    printf("Using %ld processors. \n", P);

    /* SPMD part */
    bspsieve();

    /* sequential part */
    exit(0);

} /* end main */
