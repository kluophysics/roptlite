#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int ssbtrd_(char *vect, char *uplo, integer *n, integer *kd, real *ab, integer *ldab, real *d__, real *e, real *q, integer *ldq, real *work, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper