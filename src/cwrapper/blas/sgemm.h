#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int sgemm_(char *transa, char *transb, integer *m, integer *n, integer *k, real *alpha, real *a, integer *lda, real *b, integer *ldb, real *beta, real *c__, integer *ldc);

#if __cplusplus >= 201103L
}
#endif