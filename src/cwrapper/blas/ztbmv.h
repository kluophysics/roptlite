#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int ztbmv_(char *uplo, char *trans, char *diag, integer *n, integer *k, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx);

#if __cplusplus >= 201103L
}
#endif