#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int chseqr_(char *job, char *compz, integer *n, integer *ilo, integer *ihi, complex *h__, integer *ldh, complex *w, complex *z__, integer *ldz, complex *work, integer *lwork, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper