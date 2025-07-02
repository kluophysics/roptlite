#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dtzrqf_(integer *m, integer *n, doublereal *a, integer *lda, doublereal *tau, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper