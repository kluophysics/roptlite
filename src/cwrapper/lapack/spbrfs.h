#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int spbrfs_(char *uplo, integer *n, integer *kd, integer *nrhs, real *ab, integer *ldab, real *afb, integer *ldafb, real *b, integer *ldb, real *x, integer *ldx, real *ferr, real *berr, real *work, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif