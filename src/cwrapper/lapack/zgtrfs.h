#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zgtrfs_(char *trans, integer *n, integer *nrhs, doublecomplex *dl, doublecomplex *d__, doublecomplex *du, doublecomplex *dlf, doublecomplex *df, doublecomplex *duf, doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif