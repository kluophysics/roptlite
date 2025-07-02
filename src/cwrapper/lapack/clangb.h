#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



E_f clangb_(char *norm, integer *n, integer *kl, integer *ku, complex *ab, integer *ldab, real *work);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper