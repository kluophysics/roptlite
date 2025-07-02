#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int clalsd_(char *uplo, integer *smlsiz, integer *n, integer *nrhs, real *d__, real *e, complex *b, integer *ldb, real *rcond, integer *rank, complex *work, real *rwork, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper