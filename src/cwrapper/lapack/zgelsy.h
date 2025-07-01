#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zgelsy_(integer *m, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, integer *jpvt, doublereal *rcond, integer *rank, doublecomplex *work, integer *lwork, doublereal *rwork, integer *info);

#if __cplusplus >= 201103L
}
#endif