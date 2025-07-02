#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int dgges_(char *jobvsl, char *jobvsr, char *sort, L_fp delctg, integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *sdim, doublereal *alphar, doublereal *alphai, doublereal *beta, doublereal *vsl, integer *ldvsl, doublereal *vsr, integer *ldvsr, doublereal *work, integer *lwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper