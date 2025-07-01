#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zhegv_(integer *itype, char *jobz, char *uplo, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublereal *w, doublecomplex *work, integer *lwork, doublereal *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif