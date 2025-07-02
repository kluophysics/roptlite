#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int chpgvx_(integer *itype, char *jobz, char *range, char *uplo, integer *n, complex *ap, complex *bp, real *vl, real *vu, integer *il, integer *iu, real *abstol, integer *m, real *w, complex *z__, integer *ldz, complex *work, real *rwork, integer *iwork, integer *ifail, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper