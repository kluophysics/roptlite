#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int sgges_(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, integer *n, real *a, integer *lda, real *b, integer *ldb, integer *sdim, real *alphar, real *alphai, real *beta, real *vsl, integer *ldvsl, real *vsr, integer *ldvsr, real *work, integer *lwork, logical *bwork, integer *info);

#if __cplusplus >= 201103L
}
#endif