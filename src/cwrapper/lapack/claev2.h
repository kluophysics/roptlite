#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int claev2_(complex *a, complex *b, complex *c__, real *rt1, real *rt2, real *cs1, complex *sn1);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper