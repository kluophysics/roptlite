#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dtbtrs_(char *uplo, char *trans, char *diag, integer *n, integer *kd, integer *nrhs, doublereal *ab, integer *ldab, doublereal *b, integer *ldb, integer *info);

#if __cplusplus >= 201103L
}
#endif