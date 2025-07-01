#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int cgbcon_(char *norm, integer *n, integer *kl, integer *ku, complex *ab, integer *ldab, integer *ipiv, real *anorm, real *rcond, complex *work, real *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif