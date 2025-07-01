#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dgegs_(char *jobvsl, char *jobvsr, integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *beta, doublereal *vsl, integer *ldvsl, doublereal *vsr, integer *ldvsr, doublereal *work, integer *lwork, integer *info);

#if __cplusplus >= 201103L
}
#endif