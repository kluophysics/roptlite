#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int sptsvx_(char *fact, integer *n, integer *nrhs, real *d__, real *e, real *df, real *ef, real *b, integer *ldb, real *x, integer *ldx, real *rcond, real *ferr, real *berr, real *work, integer *info);

#if __cplusplus >= 201103L
}
#endif