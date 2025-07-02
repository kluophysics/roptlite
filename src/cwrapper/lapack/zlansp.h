#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



doublereal zlansp_(char *norm, char *uplo, integer *n, doublecomplex *ap, doublereal *work);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper