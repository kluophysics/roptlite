#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int clacgv_(integer *n, complex *x, integer *incx);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper