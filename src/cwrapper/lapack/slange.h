#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



E_f slange_(char *norm, integer *m, integer *n, real *a, integer *lda, real *work);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper