#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int clarfb_(char *side, char *trans, char *direct, char *storev, integer *m, integer *n, integer *k, complex *v, integer *ldv, complex *t, integer *ldt, complex *c__, integer *ldc, complex *work, integer *ldwork);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper