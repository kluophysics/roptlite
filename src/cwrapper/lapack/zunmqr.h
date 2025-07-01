#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int zunmqr_(char *side, char *trans, integer *m, integer *n, integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork, integer *info);

#if __cplusplus >= 201103L
}
#endif