#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int slatrz_(integer *m, integer *n, integer *l, real *a, integer *lda, real *tau, real *work);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper