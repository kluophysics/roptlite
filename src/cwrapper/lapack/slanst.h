#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



E_f slanst_(char *norm, integer *n, real *d__, real *e);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper