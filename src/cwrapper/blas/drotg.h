#pragma once

#include "f2c.h" 

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  




int drotg_(doublereal *da, doublereal *db, doublereal *c__, doublereal *s);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper