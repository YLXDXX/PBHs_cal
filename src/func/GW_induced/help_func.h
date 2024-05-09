#ifndef __PBHS_GW_INDUCED_HELP_FUNC__   /* Include guard */
#define __PBHS_GW_INDUCED_HELP_FUNC__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../general.h"
#include "../other.h"

//与功率谱无关的积分函数 I(v,u,x)

int GW_I_c_func(arb_t res, const arb_t x, const arb_t y, const arb_t k, slong prec); //数值解
int GW_I_s_func(arb_t res, const arb_t x, const arb_t y, const arb_t k, slong prec);

int GW_I_c_func_analyze(arb_t res, const arb_t x, const arb_t y, const arb_t k, slong prec); //解析解
int GW_I_s_func_analyze(arb_t res, const arb_t x, const arb_t y, const arb_t k, slong prec);


int I_func_v_u_x(arb_t res, const arb_t v, const arb_t u, const arb_t x, slong prec);


#endif // __PBHS_GW_INDUCED_HELP_FUNC__  


