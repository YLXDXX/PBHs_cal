#ifndef __PBHS_PS_METHOD_PS_ABUNDANCE_ALL_K__   /* Include guard */
#define __PBHS_PS_METHOD_PS_ABUNDANCE_ALL_K__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../general.h"

#include "covariance.h"
#include "zeta_probability.h"
#include "C_l_probability.h"
#include "ps_generate_mass.h"
#include "ps_abundance.h"


int PS_abundance_beta_delta_k(arb_t res, const arb_t k_0, slong prec);

int PS_abundance_beta_delta_k_all(arb_t res, slong prec);

int PS_abundance_beta_delta_k_M(arb_t res, const arb_t k_M, slong prec);

#endif // __PBHS_PS_METHOD_PS_ABUNDANCE_ALL_K__  

