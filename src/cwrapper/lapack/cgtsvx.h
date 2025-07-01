#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int cgtsvx_(char *fact, char *trans, integer *n, integer *nrhs, complex *dl, complex *d__, complex *du, complex *dlf, complex *df, complex *duf, complex *du2, integer *ipiv, complex *b, integer *ldb, complex *x, integer *ldx, real *rcond, real *ferr, real *berr, complex *work, real *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif