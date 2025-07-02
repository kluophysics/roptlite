#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




C_f cdotc_(complex * ret_val, integer *n, complex *cx, integer *incx, complex *cy, integer *incy);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper