#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int cppsvx_(char *fact, char *uplo, integer *n, integer *nrhs, complex *ap, complex *afp, char *equed, real *s, complex *b, integer *ldb, complex *x, integer *ldx, real *rcond, real *ferr, real *berr, complex *work, real *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif