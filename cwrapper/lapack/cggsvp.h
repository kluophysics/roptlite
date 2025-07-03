#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

void roptlite_cggsvp_(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, complex *a, integer *lda, complex *b, integer *ldb, real *tola, real *tolb, integer *k, integer *l, complex *u, integer *ldu, complex *v, integer *ldv, complex *q, integer *ldq, integer *iwork, real *rwork, complex *tau, complex *work, integer *info);

#ifdef __cplusplus
}
#endif