#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dpprfs_(char *uplo, integer *n, integer *nrhs, doublereal *ap, doublereal *afp, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif