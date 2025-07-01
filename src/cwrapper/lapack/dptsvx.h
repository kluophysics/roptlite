#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dptsvx_(char *fact, integer *n, integer *nrhs, doublereal *d__, doublereal *e, doublereal *df, doublereal *ef, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *info);

#if __cplusplus >= 201103L
}
#endif