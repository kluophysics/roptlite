#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int zrotg_(doublecomplex *ca, doublecomplex *cb, doublereal *c__, doublecomplex *s);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper