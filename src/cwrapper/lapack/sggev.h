#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int sggev_(char *jobvl, char *jobvr, integer *n, real *a, integer *lda, real *b, integer *ldb, real *alphar, real *alphai, real *beta, real *vl, integer *ldvl, real *vr, integer *ldvr, real *work, integer *lwork, integer *info);

#if __cplusplus >= 201103L
}
#endif