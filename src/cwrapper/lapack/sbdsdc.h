#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int sbdsdc_(char *uplo, char *compq, integer *n, real *d__, real *e, real *u, integer *ldu, real *vt, integer *ldvt, real *q, integer *iq, real *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper