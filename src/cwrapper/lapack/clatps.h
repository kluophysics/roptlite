#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int clatps_(char *uplo, char *trans, char *diag, char *normin, integer *n, complex *ap, complex *x, real *scale, real *cnorm, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper