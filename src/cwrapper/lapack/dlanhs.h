#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



doublereal dlanhs_(char *norm, integer *n, doublereal *a, integer *lda, doublereal *work);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper