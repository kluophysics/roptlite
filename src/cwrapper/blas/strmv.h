#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int strmv_(char *uplo, char *trans, char *diag, integer *n, real *a, integer *lda, real *x, integer *incx);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper