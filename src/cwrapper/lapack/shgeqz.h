#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int shgeqz_(char *job, char *compq, char *compz, integer *n, integer *ilo, integer *ihi, real *a, integer *lda, real *b, integer *ldb, real *alphar, real *alphai, real *beta, real *q, integer *ldq, real *z__, integer *ldz, real *work, integer *lwork, integer *info);

#if __cplusplus >= 201103L
}
#endif