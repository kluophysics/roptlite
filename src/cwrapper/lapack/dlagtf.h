#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dlagtf_(integer *n, doublereal *a, doublereal *lambda, doublereal *b, doublereal *c__, doublereal *tol, doublereal *d__, integer *in, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper