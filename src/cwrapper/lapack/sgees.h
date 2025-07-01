#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int sgees_(char *jobvs, char *sort, L_fp select, integer *n, real *a, integer *lda, integer *sdim, real *wr, real *wi, real *vs, integer *ldvs, real *work, integer *lwork, logical *bwork, integer *info);

#if __cplusplus >= 201103L
}
#endif