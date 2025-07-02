#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int slasd1_(integer *nl, integer *nr, integer *sqre, real *d__, real *alpha, real *beta, real *u, integer *ldu, real *vt, integer *ldvt, integer *idxq, integer *iwork, real *work, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper