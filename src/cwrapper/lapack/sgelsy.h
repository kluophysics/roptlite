#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int sgelsy_(integer *m, integer *n, integer *nrhs, real *a, integer *lda, real *b, integer *ldb, integer *jpvt, real *rcond, integer *rank, real *work, integer *lwork, integer *info);

#if __cplusplus >= 201103L
}
#endif