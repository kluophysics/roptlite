#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int zcopy_(integer *n, doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper