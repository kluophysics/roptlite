#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dpptrf_(char *uplo, integer *n, doublereal *ap, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper