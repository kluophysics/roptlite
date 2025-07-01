#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dgeesx_(char *jobvs, char *sort, L_fp select, char *sense, integer *n, doublereal *a, integer *lda, integer *sdim, doublereal *wr, doublereal *wi, doublereal *vs, integer *ldvs, doublereal *rconde, doublereal *rcondv, doublereal *work, integer *lwork, integer *iwork, integer *liwork, logical *bwork, integer *info);

#if __cplusplus >= 201103L
}
#endif