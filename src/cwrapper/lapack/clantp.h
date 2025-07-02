#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



E_f clantp_(char *norm, char *uplo, char *diag, integer *n, complex *ap, real *work);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper