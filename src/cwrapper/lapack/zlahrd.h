#pragma once

#include "f2c.h"

namespace cwrapper {
#ifdef __cplusplus
extern "C" { 
#endif  



int zlahrd_(integer *n, integer *k, integer *nb, doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *t, integer *ldt, doublecomplex *y, integer *ldy);

#ifdef __cplusplus
}
#endif
} // end of namespace cwrapper