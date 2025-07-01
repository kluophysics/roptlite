#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int cgelsd_(integer *m, integer *n, integer *nrhs, complex *a, integer *lda, complex *b, integer *ldb, real *s, real *rcond, integer *rank, complex *work, integer *lwork, real *rwork, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif