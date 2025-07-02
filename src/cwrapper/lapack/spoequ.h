#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int spoequ_(integer *n, real *a, integer *lda, real *s, real *scond, real *amax, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper