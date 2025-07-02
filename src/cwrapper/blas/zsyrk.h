#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int zsyrk_(char *uplo, char *trans, integer *n, integer *k, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *beta, doublecomplex *c__, integer *ldc);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper