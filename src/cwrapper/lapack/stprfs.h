#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int stprfs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, real *ap, real *b, integer *ldb, real *x, integer *ldx, real *ferr, real *berr, real *work, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif