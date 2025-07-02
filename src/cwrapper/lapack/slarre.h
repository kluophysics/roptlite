#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int slarre_(integer *n, real *d__, real *e, real *tol, integer *nsplit, integer *isplit, integer *m, real *w, real *woff, real *gersch, real *work, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper