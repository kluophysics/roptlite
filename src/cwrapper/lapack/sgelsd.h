#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int sgelsd_(integer *m, integer *n, integer *nrhs, real *a, integer *lda, real *b, integer *ldb, real *s, real *rcond, integer *rank, real *work, integer *lwork, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif