#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int zlacon_(integer *n, doublecomplex *v, doublecomplex *x, doublereal *est, integer *kase);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper