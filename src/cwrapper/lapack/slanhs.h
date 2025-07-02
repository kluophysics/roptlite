#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



E_f slanhs_(char *norm, integer *n, real *a, integer *lda, real *work);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper