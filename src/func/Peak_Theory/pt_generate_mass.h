#ifndef __PBHS_PEAK_THEORY_GENERATE_MASS__   /* Include guard */
#define __PBHS_PEAK_THEORY_GENERATE_MASS__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../general.h"


int Horizon_reentry_k_to_M_H(arb_t res,const arb_t k,slong prec); //在视界重进入时，对应的视界质量为 M_H(k)

int Horizon_reentry_mu_to_M(arb_t res, const arb_t mu, const arb_t zeta_k, slong prec); //在视界重进入时，形成的黑洞的质量为 M(µ)
int Horizon_reentry_mu_to_M_relative(arb_t res, const arb_t mu, const arb_t zeta_k, slong prec); //在视界重进入时，视界质量 M_H，形成的黑洞的质量为 M，两者之比

int Horizon_reentry_M_to_mu(arb_t res,const arb_t M, const arb_t zeta_k, slong prec); //在视界重进入时，形成的黑洞的质量为M(µ)的反函数 µ(M)

int Horizon_reentry_derivative_ln_M_mu(arb_t res, const arb_t mu, const arb_t zeta_k, slong prec); //质量 M(mu) 的导数 dln(M)/dµ



#endif // __PBHS_PEAK_THEORY_GENERATE_MASS__ 
