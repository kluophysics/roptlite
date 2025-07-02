#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int sladiv_(real *a, real *b, real *c__, real *d__, real *p, real *q);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper