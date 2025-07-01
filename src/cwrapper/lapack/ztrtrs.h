#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int ztrtrs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, integer *info);

#if __cplusplus >= 201103L
}
#endif