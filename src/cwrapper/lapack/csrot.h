#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int csrot_(integer *n, complex *cx, integer *incx, complex *cy, integer *incy, real *c__, real *s);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper