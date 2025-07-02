#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



E_f slangt_(char *norm, integer *n, real *dl, real *d__, real *du);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper