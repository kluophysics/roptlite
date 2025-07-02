#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int slalsd_(char *uplo, integer *smlsiz, integer *n, integer *nrhs, real *d__, real *e, real *b, integer *ldb, real *rcond, integer *rank, real *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper