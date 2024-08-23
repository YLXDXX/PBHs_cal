#ifndef __PBHS_PS_METHOD_PS_ABUNDANCE_SIMPLE__   /* Include guard */
#define __PBHS_PS_METHOD_PS_ABUNDANCE_SIMPLE__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../general.h"

#include "covariance.h"
#include "zeta_probability.h"
#include "C_l_probability.h"
#include "ps_generate_mass.h"

struct Window_func_and_R
{
    //利用R和窗口函数求方差时用
    arb_t R;
    enum WINDOW_FUNC_TYPE w_type;
    enum WINDOW_FUNC_TYPE w_type_second;
    int method;
};


void PS_variance_help_M_to_R(arb_t res, const arb_t m, slong prec); //从PBHs的质量 M 得到相应的视界尺度 R 
void power_spectrum_density_contrast(arb_t res, const arb_t k, const arb_t R, slong prec); //密度扰动的功率谱，以 ln(k) 作为变量
void PS_variance_with_window_func(arb_t res, const arb_t R,
                                  const enum WINDOW_FUNC_TYPE w_type, const enum WINDOW_FUNC_TYPE w_type_second,
                                  const int method, slong prec); //利用窗口函数，得到质量 M 对应的统计量方差


int PS_abundance_simpele_zeta_beta_m(arb_t res, const arb_t m, slong prec); //利用曲率扰动 ζ 估算，高斯情况
int PS_abundance_simpele_zeta_beta_all(arb_t res, slong prec);
int PS_abundance_simpele_zeta_f_m(arb_t res, const arb_t m, slong prec);
int PS_abundance_simpele_zeta_f_all(arb_t res, slong prec);

int PS_abundance_simpele_delta_beta_m(arb_t res, const arb_t m, slong prec); //利用密度扰动 δ 估算，高斯情况
int PS_abundance_simpele_delta_beta_all(arb_t res, slong prec);
int PS_abundance_simpele_delta_f_m(arb_t res, const arb_t m, slong prec);
int PS_abundance_simpele_delta_f_all(arb_t res, slong prec);


int PS_abundance_simpele_compact_beta_m(arb_t res, const arb_t m, slong prec); //利用compaction function C 估算，高斯/非高斯
int PS_abundance_simpele_compact_beta_all(arb_t res, slong prec);
int PS_abundance_simpele_compact_f_m(arb_t res, const arb_t m, slong prec);
int PS_abundance_simpele_compact_f_all(arb_t res, slong prec);


#endif // __PBHS_PS_METHOD_PS_ABUNDANCE_SIMPLE__  
