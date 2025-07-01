#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dsbevd_(char *jobz, char *uplo, integer *n, integer *kd, doublereal *ab, integer *ldab, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

#if __cplusplus >= 201103L
}
#endif