#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int chetf2_(char *uplo, integer *n, complex *a, integer *lda, integer *ipiv, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper