#ifndef __PBHS_GENERAL_THRESHPLD__   /* Include guard */
#define __PBHS_GENERAL_THRESHPLD__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"

#include "basis.h"
#include "non_gaussianity.h"
#include "typical_profile.h"
#include "compaction_func.h"


int Q_parameter(arb_t res, const arb_t r_m, slong prec); //计算 q 参数，其为用来计算振幅的临界值

int Delta_c_q_parameter_simple(arb_t res, const arb_t q, slong prec); //利用 q 参数，计算临界值 δ_c（q）

int Delta_c_q_parameter_new(arb_t res, const arb_t q, slong prec); //利用 q 参数，计算临界值 δ_c（q）

int C_m_average(arb_t res, const arb_t r_m, slong prec); //求C_m的平均值

int Find_Mu_2_th(arb_t res, slong prec); //找Mu_2的临界值

int Trans_C_to_C_l(arb_t res, const arb_t x, slong prec); //已知阈值C_th时，求阈值C_{l,th}


#endif // __PBHS_GENERAL_THRESHPLD__ 
