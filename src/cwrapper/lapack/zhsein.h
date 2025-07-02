#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int zhsein_(char *side, char *eigsrc, char *initv, logical *select, integer *n, doublecomplex *h__, integer *ldh, doublecomplex *w, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr, integer *mm, integer *m, doublecomplex *work, doublereal *rwork, integer *ifaill, integer *ifailr, integer *info);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper