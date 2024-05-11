#ifndef __PBHS_GW_INDUCED_FIRST_HELP_FUNC__   /* Include guard */
#define __PBHS_GW_INDUCED_FIRST_HELP_FUNC__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../general.h"
#include "../other.h"

//诱导引力波中与功率谱无关的积分函数
//方法一 1804.07732
int GW_I_c_func(arb_t res, const arb_t x, const arb_t y, const arb_t k, slong prec); //数值解
int GW_I_s_func(arb_t res, const arb_t x, const arb_t y, const arb_t k, slong prec);

int GW_I_c_func_analyze(arb_t res, const arb_t x, const arb_t y, const arb_t k, slong prec); //解析解
int GW_I_s_func_analyze(arb_t res, const arb_t x, const arb_t y, const arb_t k, slong prec);


#endif // __PBHS_GW_INDUCED_FIRST_HELP_FUNC__  


