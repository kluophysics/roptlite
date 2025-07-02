#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int cggqrf_(integer *n, integer *m, integer *p, complex *a, integer *lda, complex *taua, complex *b, integer *ldb, complex *taub, complex *work, integer *lwork, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper