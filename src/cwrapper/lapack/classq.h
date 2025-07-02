#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int classq_(integer *n, complex *x, integer *incx, real *scale, real *sumsq);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper