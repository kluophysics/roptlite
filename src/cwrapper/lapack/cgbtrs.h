#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int cgbtrs_(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, complex *ab, integer *ldab, integer *ipiv, complex *b, integer *ldb, integer *info);

#if __cplusplus >= 201103L
}
#endif