#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int slatbs_(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd, real *ab, integer *ldab, real *x, real *scale, real *cnorm, integer *info);

#if __cplusplus >= 201103L
}
#endif