#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dsysv_(char *uplo, integer *n, integer *nrhs, doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *ldb, doublereal *work, integer *lwork, integer *info);

#if __cplusplus >= 201103L
}
#endif