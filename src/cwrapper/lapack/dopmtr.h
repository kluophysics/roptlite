#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dopmtr_(char *side, char *uplo, char *trans, integer *m, integer *n, doublereal *ap, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper