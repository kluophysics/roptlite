#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int sgeevx_(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, real *a, integer *lda, real *wr, real *wi, real *vl, integer *ldvl, real *vr, integer *ldvr, integer *ilo, integer *ihi, real *scale, real *abnrm, real *rconde, real *rcondv, real *work, integer *lwork, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif