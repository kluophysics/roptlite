#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int slasd5_(integer *i__, real *d__, real *z__, real *delta, real *rho, real *dsigma, real *work);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper