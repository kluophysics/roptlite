#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zgbsvx_(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, doublecomplex *ab, integer *ldab, doublecomplex *afb, integer *ldafb, integer *ipiv, char *equed, doublereal *r__, doublereal *c__, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif