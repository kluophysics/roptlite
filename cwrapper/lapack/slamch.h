#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

E_f roptlite_slamch_(char *cmach);
void roptlite_slamc1_(integer *beta, integer *t, logical *rnd, logical *ieee1);
void roptlite_slamc2_(integer *beta, integer *t, logical *rnd, real *eps, integer *emin, real *rmin, integer *emax, real *rmax);
E_f roptlite_slamc3_(real *a, real *b);
void roptlite_slamc4_(integer *emin, real *start, integer *base);
void roptlite_slamc5_(integer *beta, integer *p, integer *emin, logical *ieee, integer *emax, real *rmax);

#ifdef __cplusplus
}
#endif