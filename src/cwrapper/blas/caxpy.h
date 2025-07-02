#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int caxpy_(integer *n, complex *ca, complex *cx, integer *incx, complex *cy, integer *incy);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper