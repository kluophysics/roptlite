#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int claqhp_(char *uplo, integer *n, complex *ap, real *s, real *scond, real *amax, char *equed);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper