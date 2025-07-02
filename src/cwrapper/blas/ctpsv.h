#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int ctpsv_(char *uplo, char *trans, char *diag, integer *n, complex *ap, complex *x, integer *incx);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper