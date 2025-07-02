#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int sspr2_(char *uplo, integer *n, real *alpha, real *x, integer *incx, real *y, integer *incy, real *ap);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper