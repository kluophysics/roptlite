#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dgbmv_(char *trans, integer *m, integer *n, integer *kl, integer *ku, doublereal *alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, doublereal *beta, doublereal *y, integer *incy);

#if __cplusplus >= 201103L
}
#endif