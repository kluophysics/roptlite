#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int cpbcon_(char *uplo, integer *n, integer *kd, complex *ab, integer *ldab, real *anorm, real *rcond, complex *work, real *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif