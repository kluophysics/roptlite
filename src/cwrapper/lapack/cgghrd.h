#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int cgghrd_(char *compq, char *compz, integer *n, integer *ilo, integer *ihi, complex *a, integer *lda, complex *b, integer *ldb, complex *q, integer *ldq, complex *z__, integer *ldz, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper