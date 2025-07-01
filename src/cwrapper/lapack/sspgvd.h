#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int sspgvd_(integer *itype, char *jobz, char *uplo, integer *n, real *ap, real *bp, real *w, real *z__, integer *ldz, real *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

#if __cplusplus >= 201103L
}
#endif