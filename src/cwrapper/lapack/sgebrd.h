#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int sgebrd_(integer *m, integer *n, real *a, integer *lda, real *d__, real *e, real *tauq, real *taup, real *work, integer *lwork, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper