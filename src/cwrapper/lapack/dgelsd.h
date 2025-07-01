#if __cplusplus >= 201103L
extern "C" { 
#endif  

#include "f2c.h" 

int dgelsd_(integer *m, integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork, integer *iwork, integer *info);

#if __cplusplus >= 201103L
}
#endif