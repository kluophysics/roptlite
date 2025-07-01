#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dgeev_(char *jobvl, char *jobvr, integer *n, doublereal *a, integer *lda, doublereal *wr, doublereal *wi, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, doublereal *work, integer *lwork, integer *info);

#if __cplusplus >= 201103L
}
#endif