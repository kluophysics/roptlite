#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int cgtcon_(char *norm, integer *n, complex *dl, complex *d__, complex *du, complex *du2, integer *ipiv, real *anorm, real *rcond, complex *work, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper