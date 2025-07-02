#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int dscal_(integer *n, doublereal *da, doublereal *dx, integer *incx);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper