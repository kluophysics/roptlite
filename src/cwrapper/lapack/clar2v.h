#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int clar2v_(integer *n, complex *x, complex *y, complex *z__, integer *incx, real *c__, complex *s, integer *incc);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper