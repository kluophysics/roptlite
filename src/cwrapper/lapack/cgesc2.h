#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int cgesc2_(integer *n, complex *a, integer *lda, complex *rhs, integer *ipiv, integer *jpiv, real *scale);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper