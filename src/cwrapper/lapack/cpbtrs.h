#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int cpbtrs_(char *uplo, integer *n, integer *kd, integer *nrhs, complex *ab, integer *ldab, complex *b, integer *ldb, integer *info);

#if __cplusplus >= 201103L
}
#endif