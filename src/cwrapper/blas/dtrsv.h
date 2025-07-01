#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dtrsv_(char *uplo, char *trans, char *diag, integer *n, doublereal *a, integer *lda, doublereal *x, integer *incx);

#if __cplusplus >= 201103L
}
#endif