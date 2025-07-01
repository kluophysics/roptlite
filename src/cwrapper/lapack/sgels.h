#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int sgels_(char *trans, integer *m, integer *n, integer *nrhs, real *a, integer *lda, real *b, integer *ldb, real *work, integer *lwork, integer *info);

#if __cplusplus >= 201103L
}
#endif