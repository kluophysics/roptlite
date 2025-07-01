#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zheevd_(char *jobz, char *uplo, integer *n, doublecomplex *a, integer *lda, doublereal *w, doublecomplex *work, integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork, integer *liwork, integer *info);

#if __cplusplus >= 201103L
}
#endif