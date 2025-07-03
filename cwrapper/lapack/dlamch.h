#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

doublereal roptlite_dlamch_(char *cmach);
void roptlite_dlamc1_(integer *beta, integer *t, logical *rnd, logical *ieee1);
void roptlite_dlamc2_(integer *beta, integer *t, logical *rnd, doublereal *eps, integer *emin, doublereal *rmin, integer *emax, doublereal *rmax);
doublereal roptlite_dlamc3_(doublereal *a, doublereal *b);
void roptlite_dlamc4_(integer *emin, doublereal *start, integer *base);
void roptlite_dlamc5_(integer *beta, integer *p, integer *emin, logical *ieee, integer *emax, doublereal *rmax);

#ifdef __cplusplus
}
#endif