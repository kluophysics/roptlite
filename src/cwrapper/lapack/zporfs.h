#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zporfs_(char *uplo, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif