#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




doublereal dasum_(integer *n, doublereal *dx, integer *incx);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper