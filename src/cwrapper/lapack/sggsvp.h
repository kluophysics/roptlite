#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int sggsvp_(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, real *a, integer *lda, real *b, integer *ldb, real *tola, real *tolb, integer *k, integer *l, real *u, integer *ldu, real *v, integer *ldv, real *q, integer *ldq, integer *iwork, real *tau, real *work, integer *info);

#if __cplusplus >= 201103L
}
#endif