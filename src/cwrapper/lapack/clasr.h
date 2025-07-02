#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int clasr_(char *side, char *pivot, char *direct, integer *m, integer *n, real *c__, real *s, complex *a, integer *lda);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper