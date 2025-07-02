#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int ctptri_(char *uplo, char *diag, integer *n, complex *ap, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper