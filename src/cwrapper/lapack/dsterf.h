#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dsterf_(integer *n, doublereal *d__, doublereal *e, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper