#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dsbgvd_(char *jobz, char *uplo, integer *n, integer *ka, integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *ldbb, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

#if __cplusplus >= 201103L
}
#endif