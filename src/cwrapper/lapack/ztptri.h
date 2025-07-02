#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int ztptri_(char *uplo, char *diag, integer *n, doublecomplex *ap, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper