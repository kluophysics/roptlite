#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int slarfx_(char *side, integer *m, integer *n, real *v, real *tau, real *c__, integer *ldc, real *work);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper