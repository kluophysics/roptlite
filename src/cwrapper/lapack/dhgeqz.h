#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dhgeqz_(char *job, char *compq, char *compz, integer *n, integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *beta, doublereal *q, integer *ldq, doublereal *z__, integer *ldz, doublereal *work, integer *lwork, integer *info);

#if __cplusplus >= 201103L
}
#endif