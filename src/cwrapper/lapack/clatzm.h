#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int clatzm_(char *side, integer *m, integer *n, complex *v, integer *incv, complex *tau, complex *c1, complex *c2, integer *ldc, complex *work);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper