#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int csyrk_(char *uplo, char *trans, integer *n, integer *k, complex *alpha, complex *a, integer *lda, complex *beta, complex *c__, integer *ldc);

#if __cplusplus >= 201103L
}
#endif