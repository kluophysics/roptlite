#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zlatbs_(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd, doublecomplex *ab, integer *ldab, doublecomplex *x, doublereal *scale, doublereal *cnorm, integer *info);

#if __cplusplus >= 201103L
}
#endif