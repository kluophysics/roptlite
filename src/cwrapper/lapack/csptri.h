#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int csptri_(char *uplo, integer *n, complex *ap, integer *ipiv, complex *work, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper