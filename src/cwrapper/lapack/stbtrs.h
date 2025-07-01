#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int stbtrs_(char *uplo, char *trans, char *diag, integer *n, integer *kd, integer *nrhs, real *ab, integer *ldab, real *b, integer *ldb, integer *info);

#if __cplusplus >= 201103L
}
#endif