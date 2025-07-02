#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int sgbbrd_(char *vect, integer *m, integer *n, integer *ncc, integer *kl, integer *ku, real *ab, integer *ldab, real *d__, real *e, real *q, integer *ldq, real *pt, integer *ldpt, real *c__, integer *ldc, real *work, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper