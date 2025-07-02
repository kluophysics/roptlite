#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




real sdot_(integer *n, real *sx, integer *incx, real *sy, integer *incy);

//E_f sdot_(integer *n, real *sx, integer *incx, real *sy, integer *incy);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper