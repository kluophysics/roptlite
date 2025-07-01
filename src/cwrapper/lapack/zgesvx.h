#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zgesvx_(char *fact, char *trans, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, integer *ipiv, char *equed, doublereal *r__, doublereal *c__, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif