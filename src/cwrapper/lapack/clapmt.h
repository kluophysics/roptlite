#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int clapmt_(logical *forwrd, integer *m, integer *n, complex *x, integer *ldx, integer *k);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper