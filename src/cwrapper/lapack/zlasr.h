#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int zlasr_(char *side, char *pivot, char *direct, integer *m, integer *n, doublereal *c__, doublereal *s, doublecomplex *a, integer *lda);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper