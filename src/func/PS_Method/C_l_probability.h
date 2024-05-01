#ifndef __PBHS_PS_METHOD_C_l_PROBABILITY__   /* Include guard */
#define __PBHS_PS_METHOD_C_l_PROBABILITY__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../general.h"

#include "covariance.h"
#include "zeta_probability.h"


int Probability_C_l(arb_t res, const arb_t C_l, slong prec); //计算 C_ℓ 的概率密度分布 P(C_l)

int Probability_C(arb_t res, const arb_t C, slong prec); //计算compaction function C 的概率密度分布

int probability_gauss_2D(arb_t res, const arb_t x, const arb_t y, slong prec);



#endif // __PBHS_PS_METHOD_C_l_PROBABILITY__  
 
