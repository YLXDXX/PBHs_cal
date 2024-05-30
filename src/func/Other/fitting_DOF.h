#ifndef __PBHS_OTHER_FITTING_DOF__   /* Include guard */
#define __PBHS_OTHER_FITTING_DOF__

#include <arb.h> //高精度实数运算
#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"

extern arb_ptr Fitting_a_i,Fitting_b_i,Fitting_c_i,Fitting_d_i; //自由度数拟合系数
extern arb_t Fitting_m_e,Fitting_m_mu,Fitting_m_pi_0,Fitting_m_pi_pm,
             Fitting_m_1,Fitting_m_2,Fitting_m_3,Fitting_m_4;//自由度数拟合系数
extern arb_ptr Fitting_func_f_rho,Fitting_func_b_rho,Fitting_func_f_s,Fitting_func_b_s,Fitting_func_S_fit;//自由度数拟合函数系数
extern arb_ptr Fitting_func_g_star_r,Fitting_func_g_star_s;



//由温度 T Gev/Mev 拟合对应的自由度数
void Effective_degrees_of_freedom_fit(arb_t res_r, arb_t res_s, const arb_t T, char * unit, slong prec);

//波数 k Mpc^-1 对应的自由度数
void Func_k_to_degrees_of_freedom(arb_t res_r, arb_t res_s, const arb_t k, slong prec);


#endif // __PBHS_OTHER_FITTING_DOF__ 
