#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dlaeda_(integer *n, integer *tlvls, integer *curlvl, integer *curpbm, integer *prmptr, integer *perm, integer *givptr, integer *givcol, doublereal *givnum, doublereal *q, integer *qptr, doublereal *z__, doublereal *ztemp, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper