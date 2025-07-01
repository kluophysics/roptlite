#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zhpmv_(char *uplo, integer *n, doublecomplex *alpha, doublecomplex *ap, doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *incy);

#if __cplusplus >= 201103L
}
#endif