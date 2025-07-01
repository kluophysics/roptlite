#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int spbcon_(char *uplo, integer *n, integer *kd, real *ab, integer *ldab, real *anorm, real *rcond, real *work, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif