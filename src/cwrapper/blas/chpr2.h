#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int chpr2_(char *uplo, integer *n, complex *alpha, complex *x, integer *incx, complex *y, integer *incy, complex *ap);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper