#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int cbdsqr_(char *uplo, integer *n, integer *ncvt, integer *nru, integer *ncc, real *d__, real *e, complex *vt, integer *ldvt, complex *u, integer *ldu, complex *c__, integer *ldc, real *rwork, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper