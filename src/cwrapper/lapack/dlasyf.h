#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dlasyf_(char *uplo, integer *n, integer *nb, integer *kb, doublereal *a, integer *lda, integer *ipiv, doublereal *w, integer *ldw, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper