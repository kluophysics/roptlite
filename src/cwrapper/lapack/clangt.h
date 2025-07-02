#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



E_f clangt_(char *norm, integer *n, complex *dl, complex *d__, complex *du);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper