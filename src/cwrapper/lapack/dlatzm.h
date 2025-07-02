#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dlatzm_(char *side, integer *m, integer *n, doublereal *v, integer *incv, doublereal *tau, doublereal *c1, doublereal *c2, integer *ldc, doublereal *work);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper