#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int drscl_(integer *n, doublereal *sa, doublereal *sx, integer *incx);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper