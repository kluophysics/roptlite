#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int sscal_(integer *n, real *sa, real *sx, integer *incx);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper