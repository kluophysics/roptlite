#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int spbtrf_(char *uplo, integer *n, integer *kd, real *ab, integer *ldab, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper