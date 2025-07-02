#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int cung2r_(integer *m, integer *n, integer *k, complex *a, integer *lda, complex *tau, complex *work, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper