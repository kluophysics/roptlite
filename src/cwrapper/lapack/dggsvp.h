#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dggsvp_(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *tola, doublereal *tolb, integer *k, integer *l, doublereal *u, integer *ldu, doublereal *v, integer *ldv, doublereal *q, integer *ldq, integer *iwork, doublereal *tau, doublereal *work, integer *info);

#if __cplusplus >= 201103L
}
#endif