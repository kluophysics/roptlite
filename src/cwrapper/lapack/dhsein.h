#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dhsein_(char *side, char *eigsrc, char *initv, logical *select, integer *n, doublereal *h__, integer *ldh, doublereal *wr, doublereal *wi, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, integer *mm, integer *m, doublereal *work, integer *ifaill, integer *ifailr, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper