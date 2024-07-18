#ifndef __PBHS_PEAK_THEORY_NUMBER_DENSITY__   /* Include guard */
#define __PBHS_PEAK_THEORY_NUMBER_DENSITY__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../general.h"


#include "generate_mass.h"


int N_pk_help_f_xi(arb_t res, const arb_t xi, slong prec); //峰数密度相关函数 f(ξ)


int PBH_number_density_M(arb_t res,const arb_t M, slong prec); //计算质量为M的原初黑洞数密度



//峰的数密度，依照高斯类型来计算的 n_pk(mu,k)
//一般情问下，是mu和k的函数，需通过积分把k积掉
//在delta谱的情况下，大为化简，由于P_1^n 里出现delta函数，可以得到 n_pk(mu)的解析式
int Peak_number_density(arb_t res, const arb_t mu, const arb_t k, slong prec);





#endif // __PBHS_PEAK_THEORY_NUMBER_DENSITY__
