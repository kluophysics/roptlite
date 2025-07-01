#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int chbmv_(char *uplo, integer *n, integer *k, complex *alpha, complex *a, integer *lda, complex *x, integer *incx, complex *beta, complex *y, integer *incy);

#if __cplusplus >= 201103L
}
#endif