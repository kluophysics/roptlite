#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int sgbsvx_(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, real *ab, integer *ldab, real *afb, integer *ldafb, integer *ipiv, char *equed, real *r__, real *c__, real *b, integer *ldb, real *x, integer *ldx, real *rcond, real *ferr, real *berr, real *work, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif