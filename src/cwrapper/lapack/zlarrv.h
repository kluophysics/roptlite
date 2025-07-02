#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int zlarrv_(integer *n, doublereal *d__, doublereal *l, integer *isplit, integer *m, doublereal *w, integer *iblock, doublereal *gersch, doublereal *tol, doublecomplex *z__, integer *ldz, integer *isuppz, doublereal *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper