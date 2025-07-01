#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dormqr_(char *side, char *trans, integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work, integer *lwork, integer *info);

#if __cplusplus >= 201103L
}
#endif