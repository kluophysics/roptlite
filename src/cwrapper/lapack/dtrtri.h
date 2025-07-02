#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dtrtri_(char *uplo, char *diag, integer *n, doublereal *a, integer *lda, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper