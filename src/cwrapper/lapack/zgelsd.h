#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zgelsd_(integer *m, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublereal *s, doublereal *rcond, integer *rank, doublecomplex *work, integer *lwork, doublereal *rwork, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif