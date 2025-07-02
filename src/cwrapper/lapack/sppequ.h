#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int sppequ_(char *uplo, integer *n, real *ap, real *s, real *scond, real *amax, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper