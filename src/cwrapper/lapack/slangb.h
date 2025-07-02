#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



E_f slangb_(char *norm, integer *n, integer *kl, integer *ku, real *ab, integer *ldab, real *work);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper