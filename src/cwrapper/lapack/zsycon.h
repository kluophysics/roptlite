#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zsycon_(char *uplo, integer *n, doublecomplex *a, integer *lda, integer *ipiv, doublereal *anorm, doublereal *rcond, doublecomplex *work, integer *info);

#if __cplusplus >= 201103L
}
#endif