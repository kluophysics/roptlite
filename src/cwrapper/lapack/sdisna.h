#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int sdisna_(char *job, integer *m, integer *n, real *d__, real *sep, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper