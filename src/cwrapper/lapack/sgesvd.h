#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int sgesvd_(char *jobu, char *jobvt, integer *m, integer *n, real *a, integer *lda, real *s, real *u, integer *ldu, real *vt, integer *ldvt, real *work, integer *lwork, integer *info);

#if __cplusplus >= 201103L
}
#endif