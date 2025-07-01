#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int ztbcon_(char *norm, char *uplo, char *diag, integer *n, integer *kd, doublecomplex *ab, integer *ldab, doublereal *rcond, doublecomplex *work, doublereal *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif