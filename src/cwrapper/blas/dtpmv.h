#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int dtpmv_(char *uplo, char *trans, char *diag, integer *n, doublereal *ap, doublereal *x, integer *incx);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper