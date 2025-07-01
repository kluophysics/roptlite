#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dsycon_(char *uplo, integer *n, doublereal *a, integer *lda, integer *ipiv, doublereal *anorm, doublereal *rcond, doublereal *work, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif