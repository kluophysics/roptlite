#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int sgtts2_(integer *itrans, integer *n, integer *nrhs, real *dl, real *d__, real *du, real *du2, integer *ipiv, real *b, integer *ldb);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper