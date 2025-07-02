#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dspgst_(integer *itype, char *uplo, integer *n, doublereal *ap, doublereal *bp, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper