#ifndef __PBHS_GW_INDUCED_FIRST_POWER_SPECTRA__   /* Include guard */
#define __PBHS_GW_INDUCED_FIRST_POWER_SPECTRA__

#include "../../new_type.h"
#include "../../phy_constant.h"
#include "../../math_method.h"
#include "../../general.h"
#include "../../other.h"

//方法一 1804.07732 中与诱导引力波中与功率谱无关的积分函数
int GW_I_c_func(arb_t res, const arb_t x, const arb_t y, const arb_t k, slong prec); //数值解
int GW_I_s_func(arb_t res, const arb_t x, const arb_t y, const arb_t k, slong prec);

int GW_I_c_func_analyze(arb_t res, const arb_t x, const arb_t y, const arb_t k, slong prec); //解析解
int GW_I_s_func_analyze(arb_t res, const arb_t x, const arb_t y, const arb_t k, slong prec);

//方法一 1804.07732
int GW_power_spectra_01(arb_t res, const arb_t eta, const arb_t k, slong prec); //诱导引力波的功率谱
int GW_current_energy_density_01(arb_t res, const arb_t k, slong prec); //当前的引力波能量密度



#endif // __PBHS_GW_INDUCED_FIRST_POWER_SPECTRA__  


