#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int sspgst_(integer *itype, char *uplo, integer *n, real *ap, real *bp, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper