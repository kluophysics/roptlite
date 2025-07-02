#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dlartv_(integer *n, doublereal *x, integer *incx, doublereal *y, integer *incy, doublereal *c__, doublereal *s, integer *incc);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper