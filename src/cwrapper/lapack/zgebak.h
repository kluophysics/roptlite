#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int zgebak_(char *job, char *side, integer *n, integer *ilo, integer *ihi, doublereal *scale, integer *m, doublecomplex *v, integer *ldv, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper