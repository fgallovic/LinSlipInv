#ifndef HEADERS_H
#define HEADERS_H

//Comment out for single precision
#define USE_DOUBLE

#ifdef USE_DOUBLE
#define REAL double
#else
#define REAL float
#endif

//NNLS Constants
#define MAX_ITER_LS(m, n) (((m)+(n))*2)
#define MAX_ITER_NNLS(m, n)  MAX_ITER_LS(m, n)

//Utility macros
#define MIN(a, b)((a)>(b)?(b):(a))
#define MAX(a, b)((a)>(b)?(a):(b))
#define SIGN(a)((a)>0?1:-1)

//Index macros
#define M2(i, j, jl) ((i) * (jl) + (j))
#define M3(i, j, k, jl, kl) ((i) * ((jl) * (kl)) + (j) * (kl) + (k) )

#endif
