#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dgttrs_(char *trans, integer *n, integer *nrhs, doublereal *dl, doublereal *d__, doublereal *du, doublereal *du2, integer *ipiv, doublereal *b, integer *ldb, integer *info);

#if __cplusplus >= 201103L
}
#endif