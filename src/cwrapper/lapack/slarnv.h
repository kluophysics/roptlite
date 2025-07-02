#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int slarnv_(integer *idist, integer *iseed, integer *n, real *x);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper