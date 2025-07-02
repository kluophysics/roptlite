#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int slargv_(integer *n, real *x, integer *incx, real *y, integer *incy, real *c__, integer *incc);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper