#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int dswap_(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper