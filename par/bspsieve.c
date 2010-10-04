#include "bspedupack.h"
#include <limits.h>

/* This program computes all prime numbers <= a given integer N.  
*/

// prototypes: 
int blockSize(int, int, int);
int blockLow(int, int, int);
int blockHigh(int, int, int);
int blockOwner(int, int, int);
int globalIdx(int, int, int, int);
int localIdx(int, int, int, int);
double bspip(int p, int s, int n, double *x, double *y);

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
/*
double bspip(int p, int s, int n, double *x, double *y){
    // Compute inner product of vectors x and y of length n>=0 

    double inprod, *Inprod, alpha;
    int i, t;
  
    Inprod= vecallocd(p); bsp_push_reg(Inprod,p*SZDBL);
    bsp_sync();

    inprod= 0.0;
    for (i=0; i<blockSize(p,s,n); i++){
        inprod += x[i]*y[i];
    }
    for (t=0; t<p; t++){
        bsp_put(t,&inprod,Inprod,s*SZDBL,SZDBL);
    }
    bsp_sync();

    alpha= 0.0;
    for (t=0; t<p; t++){
        alpha += Inprod[t];
    }
    bsp_pop_reg(Inprod); vecfreed(Inprod);

    return alpha;

} // end bspip 

*/
void bspmarkmultiples(int p, int s, int n, int k, int *x)
{
    // mark all multiples of k as non-prime in x

    int i;
    for (i=k; i * k  <= blockHigh(p,s,n); i++)
    {
        x[localIdx(p,s,n,i*k)] = 0; //not a prime
    }


} /* end bspmarkmultiples */

int nextPrime(int p, int s, int n, int k, int *x)
{
    // find minimal i s.t. i > k and i unmarked

    int newK = k+1;
    int local = localIdx(p,s,n,newK);
    while(local < blockHigh(p,s,n) && x[local] == 0)
        local++;

    if(local > blockHigh(p,s,n))
    {
        printf("help? no primes for proc %d?\n", s);
        return INT_MAX;
    }

    // if we get here we assume we found a prime!

    printf("next prime found: %d (local index = %d)\n", globalIdx(p,s,n,local), local);

    return globalIdx(p,s,n,local);

} /* end nextPrime */

void bspsieve(){
    
    double alpha, time0, time1;
    int *x;
    int *ks; //place for proc0 to store intermediate k's
    int p, s, n, nl, i, iglob;
    int k;   // the current largest sure-prime

    n = N+1; // copy global N and increase by 1.
             // this is so the maximum array idx == N
    
    bsp_begin(P);
    p= bsp_nprocs(); /* p = number of processors obtained */ 
    s= bsp_pid();    /* s = processor number */ 
    if (s==0){
        if(n<0)
            bsp_abort("Error in input: n is negative");
        ks = vecalloci(p);
    }
//    bsp_push_reg(&n,SZINT);
//    bsp_sync();
//
//    bsp_get(0,&n,0,&n,SZINT); //everyone reads N from proc 0
//    bsp_sync();
//    bsp_pop_reg(&n);
    printf("proc %d thinks N (+1) = %d\n", s, n);

    nl= blockSize(p,s,n); // how big must s's block be?
    x= vecalloci(nl);
    for (i=0; i<nl; i++){
        // start by assuming everything is prime
        iglob= globalIdx(p,s,n,i);
        x[i]= iglob;
    }
    bsp_sync(); 
    time0=bsp_time();
    k = 2;
    // begin work

    while(k < sqrt(n))
    {
        bspmarkmultiples(p,s,n,k,x);
        k = nextPrime(p,s,n,k,x);

        bsp_push_reg(&k, SZINT);
        bsp_sync();

        if(s==0)
        {
            for(i=0;i<p; i++)
            {
                bsp_get(i, &k, 0, &ks[i], SZINT);
            }
        }

        bsp_sync();

        if(s==0)
        {
            k = findMinimum(p,ks);
            // broadcast the new minimum and continue processing
            for(i = 0; i< p ; i++)
                bsp_put(i, &k, &k, 0, SZINT);
        }
        bsp_pop_reg(&k);
        bsp_sync();
    }

    // end work
    bsp_sync();  
    time1=bsp_time();

    printf("Processor %d primes: \n", s); 
    for(i = 0; i < blockSize(p,s,n); i++)
        if( x[i] != 0)
            printf("  %d is prime\n", globalIdx(p,s,n,i));
//        printf(" x[%d] = %d\n", i, x[i]);

    fflush(stdout);
    if (s==0){
        printf("This took only %.6lf seconds.\n", time1-time0);
        fflush(stdout);
        vecfreei(ks);
    }

    vecfreei(x);
    bsp_end();

} /* end bspinprod */

int main(int argc, char **argv){

    bsp_init(bspsieve, argc, argv);

    /* sequential part */
    sscanf(argv[1], "%d", &N);
    printf("max prime requested = %d\n", N);
    P = bsp_nprocs(); // maximum amount of procs

    /* SPMD part */
    bspsieve();

    /* sequential part */
    exit(0);

} /* end main */
