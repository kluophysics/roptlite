#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int dsyr_(char *uplo, integer *n, doublereal *alpha, doublereal *x, integer *incx, doublereal *a, integer *lda);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper