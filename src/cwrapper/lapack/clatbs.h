#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int clatbs_(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd, complex *ab, integer *ldab, complex *x, real *scale, real *cnorm, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper