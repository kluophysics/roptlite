#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int slarz_(char *side, integer *m, integer *n, integer *l, real *v, integer *incv, real *tau, real *c__, integer *ldc, real *work);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper