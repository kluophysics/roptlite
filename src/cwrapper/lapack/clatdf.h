#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int clatdf_(integer *ijob, integer *n, complex *z__, integer *ldz, complex *rhs, real *rdsum, real *rdscal, integer *ipiv, integer *jpiv);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper