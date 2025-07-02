#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int cpbtrs_(char *uplo, integer *n, integer *kd, integer *nrhs, complex *ab, integer *ldab, complex *b, integer *ldb, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper