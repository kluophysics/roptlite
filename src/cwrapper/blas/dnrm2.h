#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




doublereal dnrm2_(integer *n, doublereal *x, integer *incx);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper