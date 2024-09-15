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

int Func_background_coupled_odes(arb_ptr yp, const arb_t x, const arb_ptr y, const slong dim,
                      void* param, const slong order, slong prec);

int Func_perturbation_phi_odes(arb_ptr yp, const arb_t t, const arb_ptr y, const slong dim,
                               void* param, const slong order, slong prec);

#endif // __PBHS_POWER_SPECTRA_NUMERIC_INFLATION_ODES_SLOVE__ 
