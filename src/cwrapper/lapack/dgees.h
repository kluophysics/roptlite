#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dgees_(char *jobvs, char *sort, L_fp select, integer *n, doublereal *a, integer *lda, integer *sdim, doublereal *wr, doublereal *wi, doublereal *vs, integer *ldvs, doublereal *work, integer *lwork, logical *bwork, integer *info);

#if __cplusplus >= 201103L
}
#endif