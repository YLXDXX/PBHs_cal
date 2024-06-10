#ifndef __PBHS_GENERAL_BASIS__   /* Include guard */
#define __PBHS_GENERAL_BASIS__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../other.h"

#include <arb.h> //高精度实数运算

void Upward_step_power_spectrum_k_c_to_k_s(arb_t k_s, const arb_t k_c, slong prec);

int window_function(arb_t res, const arb_t k, slong prec); //在k空间的窗口函数

int Linear_transfer_function(arb_t res, const arb_t k, const arb_t eta, slong prec);//线性转移函数

int power_spectrum(arb_t res, const arb_t k, slong prec); //功率谱

#endif // __PBHS_GENERAL_BASIS__ 
