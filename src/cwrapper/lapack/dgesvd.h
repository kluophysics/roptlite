#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n, doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, integer *info);

#if __cplusplus >= 201103L
}
#endif