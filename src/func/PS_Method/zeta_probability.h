#ifndef __PBHS_PS_METHOD_ZETA_PROBABILITY__   /* Include guard */
#define __PBHS_PS_METHOD_ZETA_PROBABILITY__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../general.h"

#include "covariance.h"

int interior_probability_gauss(arb_t res, const arb_t zeta, slong prec); //高斯扰动ζ的概率 P(ζ_G)=P(Y)

int Probability_zeta(arb_t res, const arb_t zeta,void* p, const slong order, slong prec); //计算 ζ 的概率密度 P(ζ)


#endif // __PBHS_PS_METHOD_ZETA_PROBABILITY__  
 
