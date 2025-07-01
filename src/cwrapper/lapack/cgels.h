#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int cgels_(char *trans, integer *m, integer *n, integer *nrhs, complex *a, integer *lda, complex *b, integer *ldb, complex *work, integer *lwork, integer *info);

#if __cplusplus >= 201103L
}
#endif