#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



integer izmax1_(integer *n, doublecomplex *cx, integer *incx);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper