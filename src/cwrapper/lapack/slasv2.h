#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int slasv2_(real *f, real *g, real *h__, real *ssmin, real *ssmax, real *snr, real *csr, real *snl, real *csl);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper