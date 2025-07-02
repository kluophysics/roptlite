#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int slarrb_(integer *n, real *d__, real *l, real *ld, real *lld, integer *ifirst, integer *ilast, real *sigma, real *reltol, real *w, real *wgap, real *werr, real *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper