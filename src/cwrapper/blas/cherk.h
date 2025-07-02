#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int cherk_(char *uplo, char *trans, integer *n, integer *k, real *alpha, complex *a, integer *lda, real *beta, complex *c__, integer *ldc);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper