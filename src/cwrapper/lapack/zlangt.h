#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



doublereal zlangt_(char *norm, integer *n, doublecomplex *dl, doublecomplex *d__, doublecomplex *du);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper