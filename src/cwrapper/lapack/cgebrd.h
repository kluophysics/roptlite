#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int cgebrd_(integer *m, integer *n, complex *a, integer *lda, real *d__, real *e, complex *tauq, complex *taup, complex *work, integer *lwork, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper