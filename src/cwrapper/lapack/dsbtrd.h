#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dsbtrd_(char *vect, char *uplo, integer *n, integer *kd, doublereal *ab, integer *ldab, doublereal *d__, doublereal *e, doublereal *q, integer *ldq, doublereal *work, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper