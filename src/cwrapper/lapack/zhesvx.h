#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zhesvx_(char *fact, char *uplo, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *work, integer *lwork, doublereal *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif