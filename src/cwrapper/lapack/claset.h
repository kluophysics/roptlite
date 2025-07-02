#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int claset_(char *uplo, integer *m, integer *n, complex *alpha, complex *beta, complex *a, integer *lda);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper