#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int cpptri_(char *uplo, integer *n, complex *ap, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper