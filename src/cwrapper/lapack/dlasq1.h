#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dlasq1_(integer *n, doublereal *d__, doublereal *e, doublereal *work, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper