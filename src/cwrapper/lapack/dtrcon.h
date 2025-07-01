#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dtrcon_(char *norm, char *uplo, char *diag, integer *n, doublereal *a, integer *lda, doublereal *rcond, doublereal *work, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif