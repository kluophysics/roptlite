#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int cpbsvx_(char *fact, char *uplo, integer *n, integer *kd, integer *nrhs, complex *ab, integer *ldab, complex *afb, integer *ldafb, char *equed, real *s, complex *b, integer *ldb, complex *x, integer *ldx, real *rcond, real *ferr, real *berr, complex *work, real *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif