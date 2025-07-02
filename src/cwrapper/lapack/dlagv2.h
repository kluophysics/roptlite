#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dlagv2_(doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *beta, doublereal *csl, doublereal *snl, doublereal *csr, doublereal *snr);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper