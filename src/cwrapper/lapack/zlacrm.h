#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int zlacrm_(integer *m, integer *n, doublecomplex *a, integer *lda, doublereal *b, integer *ldb, doublecomplex *c__, integer *ldc, doublereal *rwork);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper