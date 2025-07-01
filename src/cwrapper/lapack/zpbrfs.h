#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zpbrfs_(char *uplo, integer *n, integer *kd, integer *nrhs, doublecomplex *ab, integer *ldab, doublecomplex *afb, integer *ldafb, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif