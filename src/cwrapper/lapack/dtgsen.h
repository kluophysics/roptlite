#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dtgsen_(integer *ijob, logical *wantq, logical *wantz, logical *select, integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *beta, doublereal *q, integer *ldq, doublereal *z__, integer *ldz, integer *m, doublereal *pl, doublereal *pr, doublereal *dif, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

#if __cplusplus >= 201103L
}
#endif