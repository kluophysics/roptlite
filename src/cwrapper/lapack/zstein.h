#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int zstein_(integer *n, doublereal *d__, doublereal *e, integer *m, doublereal *w, integer *iblock, integer *isplit, doublecomplex *z__, integer *ldz, doublereal *work, integer *iwork, integer *ifail, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper