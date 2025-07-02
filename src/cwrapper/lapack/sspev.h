#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int sspev_(char *jobz, char *uplo, integer *n, real *ap, real *w, real *z__, integer *ldz, real *work, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper