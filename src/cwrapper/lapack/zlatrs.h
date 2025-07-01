#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zlatrs_(char *uplo, char *trans, char *diag, char *normin, integer *n, doublecomplex *a, integer *lda, doublecomplex *x, doublereal *scale, doublereal *cnorm, integer *info);

#if __cplusplus >= 201103L
}
#endif