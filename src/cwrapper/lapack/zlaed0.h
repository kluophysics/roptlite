#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int zlaed0_(integer *qsiz, integer *n, doublereal *d__, doublereal *e, doublecomplex *q, integer *ldq, doublecomplex *qstore, integer *ldqs, doublereal *rwork, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper