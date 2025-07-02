#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int zgtcon_(char *norm, integer *n, doublecomplex *dl, doublecomplex *d__, doublecomplex *du, doublecomplex *du2, integer *ipiv, doublereal *anorm, doublereal *rcond, doublecomplex *work, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper