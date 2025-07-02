#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int cgetf2_(integer *m, integer *n, complex *a, integer *lda, integer *ipiv, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper