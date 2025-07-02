#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int cupmtr_(char *side, char *uplo, char *trans, integer *m, integer *n, complex *ap, complex *tau, complex *c__, integer *ldc, complex *work, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper