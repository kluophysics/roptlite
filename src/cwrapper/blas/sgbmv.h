#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int sgbmv_(char *trans, integer *m, integer *n, integer *kl, integer *ku, real *alpha, real *a, integer *lda, real *x, integer *incx, real *beta, real *y, integer *incy);

#if __cplusplus >= 201103L
}
#endif