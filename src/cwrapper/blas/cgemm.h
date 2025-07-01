#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int cgemm_(char *transa, char *transb, integer *m, integer *n, integer *k, complex *alpha, complex *a, integer *lda, complex *b, integer *ldb, complex *beta, complex *c__, integer *ldc);

#if __cplusplus >= 201103L
}
#endif