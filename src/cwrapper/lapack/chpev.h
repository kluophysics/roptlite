#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int chpev_(char *jobz, char *uplo, integer *n, complex *ap, real *w, complex *z__, integer *ldz, complex *work, real *rwork, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper