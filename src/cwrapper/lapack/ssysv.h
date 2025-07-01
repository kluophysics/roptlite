#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int ssysv_(char *uplo, integer *n, integer *nrhs, real *a, integer *lda, integer *ipiv, real *b, integer *ldb, real *work, integer *lwork, integer *info);

#if __cplusplus >= 201103L
}
#endif