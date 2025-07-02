#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int clatrd_(char *uplo, integer *n, integer *nb, complex *a, integer *lda, real *e, complex *tau, complex *w, integer *ldw);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper