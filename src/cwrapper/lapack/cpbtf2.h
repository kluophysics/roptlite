#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int cpbtf2_(char *uplo, integer *n, integer *kd, complex *ab, integer *ldab, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper