#ifndef __PBHS_POWER_SPECTRA_NUMERIC_INFLATION_ODES_SLOVE__   /* Include guard */
#define __PBHS_POWER_SPECTRA_NUMERIC_INFLATION_ODES_SLOVE__

#include <arb.h> //高精度实数运算
#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"


//扰动求解
extern Interp_coe_t G_Interp_coe_0;
extern Interp_coe_t G_Interp_coe_1;
extern Interp_coe_t G_Interp_coe_2;
extern Interp_coe_t G_Interp_coe_3;
extern Interp_coe_t G_Interp_coe_4;
extern arb_t G_fourier_k;

void Func_Smoothing_Step_Function(arb_t res, const arb_t x, const arb_t smoothness, slong prec);
void Func_S_prime(arb_t res, const arb_t x, const arb_t smoothness, slong prec);
void Func_S_pp(arb_t res, const arb_t x, const arb_t smoothness, slong prec);
void Func_V_phi(arb_t res, const arb_t phi, slong prec);
void Func_V_phi_p(arb_t res, const arb_t phi, slong prec);
void Func_V_phi_pp(arb_t res, const arb_t phi, slong prec);


int Func_background_coupled_odes(arb_ptr yp, const arb_t x, const arb_ptr y, const slong dim,
                      void* param, const slong order, slong prec);

int Func_perturbation_phi_odes(arb_ptr yp, const arb_t t, const arb_ptr y, const slong dim,
                               void* param, const slong order, slong prec);

//得到暴胀中，模式 k 进入视界的时间及此时的 N
void Func_get_time_k_enter(arb_t t, arb_t N, const arb_t k, const arb_t t_a, const arb_t t_b, //在区间 [t_a, t_b] 内找根
                           const arb_t a_i, const ODEs_DOPRI54_dense_t d_out, // a_i 为尺度因子的初始值
                           slong prec);

#endif // __PBHS_POWER_SPECTRA_NUMERIC_INFLATION_ODES_SLOVE__ 
