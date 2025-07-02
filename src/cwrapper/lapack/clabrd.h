#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int clabrd_(integer *m, integer *n, integer *nb, complex *a, integer *lda, real *d__, real *e, complex *tauq, complex *taup, complex *x, integer *ldx, complex *y, integer *ldy);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper