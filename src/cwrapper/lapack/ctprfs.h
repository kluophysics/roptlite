#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int ctprfs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, complex *ap, complex *b, integer *ldb, complex *x, integer *ldx, real *ferr, real *berr, complex *work, real *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif