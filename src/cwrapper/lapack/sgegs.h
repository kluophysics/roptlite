#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int sgegs_(char *jobvsl, char *jobvsr, integer *n, real *a, integer *lda, real *b, integer *ldb, real *alphar, real *alphai, real *beta, real *vsl, integer *ldvsl, real *vsr, integer *ldvsr, real *work, integer *lwork, integer *info);

#if __cplusplus >= 201103L
}
#endif