#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zhpgvd_(integer *itype, char *jobz, char *uplo, integer *n, doublecomplex *ap, doublecomplex *bp, doublereal *w, doublecomplex *z__, integer *ldz, doublecomplex *work, integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork, integer *liwork, integer *info);

#if __cplusplus >= 201103L
}
#endif