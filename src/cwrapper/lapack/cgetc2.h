#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int cgetc2_(integer *n, complex *a, integer *lda, integer *ipiv, integer *jpiv, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper