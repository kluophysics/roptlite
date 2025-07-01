#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zgbmv_(char *trans, integer *m, integer *n, integer *kl, integer *ku, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *incy);

#if __cplusplus >= 201103L
}
#endif