#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int cspr_(char *uplo, integer *n, complex *alpha, complex *x, integer *incx, complex *ap);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper