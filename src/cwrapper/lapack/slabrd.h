#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int slabrd_(integer *m, integer *n, integer *nb, real *a, integer *lda, real *d__, real *e, real *tauq, real *taup, real *x, integer *ldx, real *y, integer *ldy);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper