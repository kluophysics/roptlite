#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int cgtts2_(integer *itrans, integer *n, integer *nrhs, complex *dl, complex *d__, complex *du, complex *du2, integer *ipiv, complex *b, integer *ldb);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper