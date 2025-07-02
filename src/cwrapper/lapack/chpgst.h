#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int chpgst_(integer *itype, char *uplo, integer *n, complex *ap, complex *bp, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper