#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int csysv_(char *uplo, integer *n, integer *nrhs, complex *a, integer *lda, integer *ipiv, complex *b, integer *ldb, complex *work, integer *lwork, integer *info);

#if __cplusplus >= 201103L
}
#endif