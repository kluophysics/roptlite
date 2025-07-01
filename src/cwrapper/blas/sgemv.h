#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int sgemv_(char *trans, integer *m, integer *n, real *alpha, real *a, integer *lda, real *x, integer *incx, real *beta, real *y, integer *incy);

#if __cplusplus >= 201103L
}
#endif