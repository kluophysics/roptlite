#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dsyevd_(char *jobz, char *uplo, integer *n, doublereal *a, integer *lda, doublereal *w, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

#if __cplusplus >= 201103L
}
#endif