#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int spteqr_(char *compz, integer *n, real *d__, real *e, real *z__, integer *ldz, real *work, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper