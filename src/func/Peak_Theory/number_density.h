#ifndef __PBHS_PEAK_THEORY_NUMBER_DENSITY__   /* Include guard */
#define __PBHS_PEAK_THEORY_NUMBER_DENSITY__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../general.h"

#include "generate_mass.h"


int N_pk_help_f_xi(arb_t res, const arb_t xi, slong prec); //峰数密度相关函数 f(ξ)

int Peak_number_density(arb_t res, const arb_t mu, const arb_t k, slong prec); //峰的数密度 n_pk(μ, k)

int PBH_number_density_M(arb_t res,const arb_t M, slong prec); //计算质量为M的原初黑洞数密度


#endif // __PBHS_PEAK_THEORY_NUMBER_DENSITY__
