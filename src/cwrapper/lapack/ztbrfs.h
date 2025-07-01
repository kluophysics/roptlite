#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int ztbrfs_(char *uplo, char *trans, char *diag, integer *n, integer *kd, integer *nrhs, doublecomplex *ab, integer *ldab, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif