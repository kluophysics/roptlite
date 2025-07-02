#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int cptts2_(integer *iuplo, integer *n, integer *nrhs, real *d__, complex *e, complex *b, integer *ldb);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper