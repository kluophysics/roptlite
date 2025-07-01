#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int ztrmm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb);

#if __cplusplus >= 201103L
}
#endif