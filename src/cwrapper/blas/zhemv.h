#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zhemv_(char *uplo, integer *n, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *incy);

#if __cplusplus >= 201103L
}
#endif