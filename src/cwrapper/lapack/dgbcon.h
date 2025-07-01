#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dgbcon_(char *norm, integer *n, integer *kl, integer *ku, doublereal *ab, integer *ldab, integer *ipiv, doublereal *anorm, doublereal *rcond, doublereal *work, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif