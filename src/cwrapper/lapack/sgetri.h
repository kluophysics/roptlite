#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int sgetri_(integer *n, real *a, integer *lda, integer *ipiv, real *work, integer *lwork, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper