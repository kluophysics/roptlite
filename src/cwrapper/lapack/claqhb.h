#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int claqhb_(char *uplo, integer *n, integer *kd, complex *ab, integer *ldab, real *s, real *scond, real *amax, char *equed);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper