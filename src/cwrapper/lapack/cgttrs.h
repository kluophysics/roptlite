#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int cgttrs_(char *trans, integer *n, integer *nrhs, complex *dl, complex *d__, complex *du, complex *du2, integer *ipiv, complex *b, integer *ldb, integer *info);

#if __cplusplus >= 201103L
}
#endif