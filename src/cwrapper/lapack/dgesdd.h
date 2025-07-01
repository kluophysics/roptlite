#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dgesdd_(char *jobz, integer *m, integer *n, doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif