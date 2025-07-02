#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int csscal_(integer *n, real *sa, complex *cx, integer *incx);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper