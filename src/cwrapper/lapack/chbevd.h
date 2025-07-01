#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int chbevd_(char *jobz, char *uplo, integer *n, integer *kd, complex *ab, integer *ldab, real *w, complex *z__, integer *ldz, complex *work, integer *lwork, real *rwork, integer *lrwork, integer *iwork, integer *liwork, integer *info);

#if __cplusplus >= 201103L
}
#endif