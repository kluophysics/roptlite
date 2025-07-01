#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int sgerfs_(char *trans, integer *n, integer *nrhs, real *a, integer *lda, real *af, integer *ldaf, integer *ipiv, real *b, integer *ldb, real *x, integer *ldx, real *ferr, real *berr, real *work, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif