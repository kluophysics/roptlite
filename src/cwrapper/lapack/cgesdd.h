#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int cgesdd_(char *jobz, integer *m, integer *n, complex *a, integer *lda, real *s, complex *u, integer *ldu, complex *vt, integer *ldvt, complex *work, integer *lwork, real *rwork, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif