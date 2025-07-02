#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int claein_(logical *rightv, logical *noinit, integer *n, complex *h__, integer *ldh, complex *w, complex *v, complex *b, integer *ldb, real *rwork, real *eps3, real *smlnum, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper