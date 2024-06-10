#ifndef __PBHS_PS_METHOD_PS_ABUNDANCE__   /* Include guard */
#define __PBHS_PS_METHOD_PS_ABUNDANCE__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../general.h"

#include "covariance.h"
#include "zeta_probability.h"
#include "C_l_probability.h"
#include "ps_generate_mass.h"


int PS_abundance_beta_m(arb_t res, const arb_t m, slong prec); //PBHs的质量分布 β(m)

int PS_abundance_beta_all(arb_t res, slong prec);

void beta_m_to_f_m_coefficient(arb_t res, slong prec); //生成时的能量密度 β(m)/β 和当今的能量密度 f(m)/f 间的转换系数

int PS_abundance_f_m(arb_t res, const arb_t m, slong prec); //原初黑洞在暗物质中关于质量分布的占比 f(m)

int PS_abundance_f_all(arb_t res, slong prec);


#endif // __PBHS_PS_METHOD_PS_ABUNDANCE__  
 
