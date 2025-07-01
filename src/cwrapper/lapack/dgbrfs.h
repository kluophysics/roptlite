#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dgbrfs_(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, doublereal *ab, integer *ldab, doublereal *afb, integer *ldafb, integer *ipiv, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif