#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int ssymm_(char *side, char *uplo, integer *m, integer *n, real *alpha, real *a, integer *lda, real *b, integer *ldb, real *beta, real *c__, integer *ldc);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper