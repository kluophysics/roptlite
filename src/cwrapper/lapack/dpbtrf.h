#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dpbtrf_(char *uplo, integer *n, integer *kd, doublereal *ab, integer *ldab, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper