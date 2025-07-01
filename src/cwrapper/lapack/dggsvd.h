#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dggsvd_(char *jobu, char *jobv, char *jobq, integer *m, integer *n, integer *p, integer *k, integer *l, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *alpha, doublereal *beta, doublereal *u, integer *ldu, doublereal *v, integer *ldv, doublereal *q, integer *ldq, doublereal *work, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif