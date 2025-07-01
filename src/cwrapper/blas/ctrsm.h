#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int ctrsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, complex *alpha, complex *a, integer *lda, complex *b, integer *ldb);

#if __cplusplus >= 201103L
}
#endif