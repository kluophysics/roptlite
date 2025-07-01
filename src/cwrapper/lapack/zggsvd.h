#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zggsvd_(char *jobu, char *jobv, char *jobq, integer *m, integer *n, integer *p, integer *k, integer *l, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublereal *alpha, doublereal *beta, doublecomplex *u, integer *ldu, doublecomplex *v, integer *ldv, doublecomplex *q, integer *ldq, doublecomplex *work, doublereal *rwork, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif