#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zgeesx_(char *jobvs, char *sort, L_fp select, char *sense, integer *n, doublecomplex *a, integer *lda, integer *sdim, doublecomplex *w, doublecomplex *vs, integer *ldvs, doublereal *rconde, doublereal *rcondv, doublecomplex *work, integer *lwork, doublereal *rwork, logical *bwork, integer *info);

#if __cplusplus >= 201103L
}
#endif