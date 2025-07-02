#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int zlacgv_(integer *n, doublecomplex *x, integer *incx);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper