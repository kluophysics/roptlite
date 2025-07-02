#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dlassq_(integer *n, doublereal *x, integer *incx, doublereal *scale, doublereal *sumsq);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper