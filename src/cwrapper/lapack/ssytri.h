#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int ssytri_(char *uplo, integer *n, real *a, integer *lda, integer *ipiv, real *work, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper