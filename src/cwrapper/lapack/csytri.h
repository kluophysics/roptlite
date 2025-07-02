#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int csytri_(char *uplo, integer *n, complex *a, integer *lda, integer *ipiv, complex *work, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper