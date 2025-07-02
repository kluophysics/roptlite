#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int slatdf_(integer *ijob, integer *n, real *z__, integer *ldz, real *rhs, real *rdsum, real *rdscal, integer *ipiv, integer *jpiv);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper