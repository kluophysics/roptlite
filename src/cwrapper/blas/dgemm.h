#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dgemm_(char *transa, char *transb, integer *m, integer *n, integer *k, doublereal *alpha, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, integer *ldc);

#if __cplusplus >= 201103L
}
#endif