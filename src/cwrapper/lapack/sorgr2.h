#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int sorgr2_(integer *m, integer *n, integer *k, real *a, integer *lda, real *tau, real *work, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper