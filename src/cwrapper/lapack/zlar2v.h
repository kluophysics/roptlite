#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int zlar2v_(integer *n, doublecomplex *x, doublecomplex *y, doublecomplex *z__, integer *incx, doublereal *c__, doublecomplex *s, integer *incc);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper