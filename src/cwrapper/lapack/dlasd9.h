#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dlasd9_(integer *icompq, integer *ldu, integer *k, doublereal *d__, doublereal *z__, doublereal *vf, doublereal *vl, doublereal *difl, doublereal *difr, doublereal *dsigma, doublereal *work, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper