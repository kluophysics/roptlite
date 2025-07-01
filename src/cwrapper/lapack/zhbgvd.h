#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zhbgvd_(char *jobz, char *uplo, integer *n, integer *ka, integer *kb, doublecomplex *ab, integer *ldab, doublecomplex *bb, integer *ldbb, doublereal *w, doublecomplex *z__, integer *ldz, doublecomplex *work, integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork, integer *liwork, integer *info);

#if __cplusplus >= 201103L
}
#endif