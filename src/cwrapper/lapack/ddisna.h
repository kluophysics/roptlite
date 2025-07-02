#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int ddisna_(char *job, integer *m, integer *n, doublereal *d__, doublereal *sep, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper