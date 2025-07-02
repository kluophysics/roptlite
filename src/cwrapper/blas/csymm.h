#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int csymm_(char *side, char *uplo, integer *m, integer *n, complex *alpha, complex *a, integer *lda, complex *b, integer *ldb, complex *beta, complex *c__, integer *ldc);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper