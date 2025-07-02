#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int ssterf_(integer *n, real *d__, real *e, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper