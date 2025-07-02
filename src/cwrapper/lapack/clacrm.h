#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int clacrm_(integer *m, integer *n, complex *a, integer *lda, real *b, integer *ldb, complex *c__, integer *ldc, real *rwork);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper