#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



E_f scsum1_(integer *n, complex *cx, integer *incx);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper