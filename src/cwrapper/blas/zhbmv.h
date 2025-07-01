#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zhbmv_(char *uplo, integer *n, integer *k, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *incy);

#if __cplusplus >= 201103L
}
#endif