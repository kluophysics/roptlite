#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int ctbcon_(char *norm, char *uplo, char *diag, integer *n, integer *kd, complex *ab, integer *ldab, real *rcond, complex *work, real *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif