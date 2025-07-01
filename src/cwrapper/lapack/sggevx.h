#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int sggevx_(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, real *a, integer *lda, real *b, integer *ldb, real *alphar, real *alphai, real *beta, real *vl, integer *ldvl, real *vr, integer *ldvr, integer *ilo, integer *ihi, real *lscale, real *rscale, real *abnrm, real *bbnrm, real *rconde, real *rcondv, real *work, integer *lwork, integer *iwork, logical *bwork, integer *info);

#if __cplusplus >= 201103L
}
#endif