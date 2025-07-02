#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int cher2k_(char *uplo, char *trans, integer *n, integer *k, complex *alpha, complex *a, integer *lda, complex *b, integer *ldb, real *beta, complex *c__, integer *ldc);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper