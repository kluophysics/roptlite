#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int cgelsx_(integer *m, integer *n, integer *nrhs, complex *a, integer *lda, complex *b, integer *ldb, integer *jpvt, real *rcond, integer *rank, complex *work, real *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif