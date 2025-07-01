#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int spprfs_(char *uplo, integer *n, integer *nrhs, real *ap, real *afp, real *b, integer *ldb, real *x, integer *ldx, real *ferr, real *berr, real *work, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif