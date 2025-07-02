#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




real snrm2_(integer *n, real *x, integer *incx);

//E_f snrm2_(integer *n, real *x, integer *incx);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper