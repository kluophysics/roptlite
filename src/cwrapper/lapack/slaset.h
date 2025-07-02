#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int slaset_(char *uplo, integer *m, integer *n, real *alpha, real *beta, real *a, integer *lda);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper