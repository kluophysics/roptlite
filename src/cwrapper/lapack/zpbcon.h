#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zpbcon_(char *uplo, integer *n, integer *kd, doublecomplex *ab, integer *ldab, doublereal *anorm, doublereal *rcond, doublecomplex *work, doublereal *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif