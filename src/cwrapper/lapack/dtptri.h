#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dtptri_(char *uplo, char *diag, integer *n, doublereal *ap, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper