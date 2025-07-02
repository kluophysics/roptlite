#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int slatrd_(char *uplo, integer *n, integer *nb, real *a, integer *lda, real *e, real *tau, real *w, integer *ldw);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper