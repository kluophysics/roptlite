#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int cupgtr_(char *uplo, integer *n, complex *ap, complex *tau, complex *q, integer *ldq, complex *work, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper