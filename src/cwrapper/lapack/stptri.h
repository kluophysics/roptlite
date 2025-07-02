#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int stptri_(char *uplo, char *diag, integer *n, real *ap, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper