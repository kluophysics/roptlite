#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



doublereal zlanht_(char *norm, integer *n, doublereal *d__, doublecomplex *e);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper