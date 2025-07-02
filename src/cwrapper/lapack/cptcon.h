#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int cptcon_(integer *n, real *d__, complex *e, real *anorm, real *rcond, real *rwork, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper