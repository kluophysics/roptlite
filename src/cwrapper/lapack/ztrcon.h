#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int ztrcon_(char *norm, char *uplo, char *diag, integer *n, doublecomplex *a, integer *lda, doublereal *rcond, doublecomplex *work, doublereal *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif