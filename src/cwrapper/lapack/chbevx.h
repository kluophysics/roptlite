#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int chbevx_(char *jobz, char *range, char *uplo, integer *n, integer *kd, complex *ab, integer *ldab, complex *q, integer *ldq, real *vl, real *vu, integer *il, integer *iu, real *abstol, integer *m, real *w, complex *z__, integer *ldz, complex *work, real *rwork, integer *iwork, integer *ifail, integer *info);

#if __cplusplus >= 201103L
}
#endif