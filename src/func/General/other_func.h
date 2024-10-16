#ifndef __PBHS_GENERAL_OTHER_FUNC__   /* Include guard */
#define __PBHS_GENERAL_OTHER_FUNC__ 

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"

#include "basis.h"
#include "non_gaussianity.h"
#include "typical_profile.h"
#include "compaction_func.h"
#include "threshold.h"

int Trans_C_to_C_l(arb_t res, const arb_t x, slong prec); //已知阈值C_th时，求阈值C_{l,th}

int Get_PK_mu_max(arb_t res, const arb_t zeta_k, slong prec); //求出参数 μ 的上限

void beta_m_to_f_m_coefficient(arb_t res, slong prec); //生成时的能量密度 β(m)/β 和当今的能量密度 f(m)/f 间的转换系数

void Get_all_k_over_k_ch(arb_t k_ch_times_r_m, arb_t k_ch, const arb_t x_m, slong prec); //考虑所有 k 模式时，特征模式的求解

void Power_spectra_convolution(arb_t res, const arb_t k, slong prec); //功率谱的卷积

#endif // __PBHS_GENERAL_OTHER_FUNC__ 
