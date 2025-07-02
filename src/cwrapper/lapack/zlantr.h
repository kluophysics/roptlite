#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



doublereal zlantr_(char *norm, char *uplo, char *diag, integer *m, integer *n, doublecomplex *a, integer *lda, doublereal *work);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper