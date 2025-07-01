#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int cunmqr_(char *side, char *trans, integer *m, integer *n, integer *k, complex *a, integer *lda, complex *tau, complex *c__, integer *ldc, complex *work, integer *lwork, integer *info);

#if __cplusplus >= 201103L
}
#endif