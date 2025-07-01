#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int chpgvd_(integer *itype, char *jobz, char *uplo, integer *n, complex *ap, complex *bp, real *w, complex *z__, integer *ldz, complex *work, integer *lwork, real *rwork, integer *lrwork, integer *iwork, integer *liwork, integer *info);

#if __cplusplus >= 201103L
}
#endif