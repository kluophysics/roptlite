#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int sgbcon_(char *norm, integer *n, integer *kl, integer *ku, real *ab, integer *ldab, integer *ipiv, real *anorm, real *rcond, real *work, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif