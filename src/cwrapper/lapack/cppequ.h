#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int cppequ_(char *uplo, integer *n, complex *ap, real *s, real *scond, real *amax, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper