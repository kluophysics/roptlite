#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int ssyr_(char *uplo, integer *n, real *alpha, real *x, integer *incx, real *a, integer *lda);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper