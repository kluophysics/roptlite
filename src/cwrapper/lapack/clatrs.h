#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int clatrs_(char *uplo, char *trans, char *diag, char *normin, integer *n, complex *a, integer *lda, complex *x, real *scale, real *cnorm, integer *info);

#if __cplusplus >= 201103L
}
#endif