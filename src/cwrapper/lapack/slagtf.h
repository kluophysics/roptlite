#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int slagtf_(integer *n, real *a, real *lambda, real *b, real *c__, real *tol, real *d__, integer *in, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper