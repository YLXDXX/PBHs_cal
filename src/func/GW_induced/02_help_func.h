#ifndef __PBHS_GW_INDUCED_SECOND_HELP_FUNC__   /* Include guard */
#define __PBHS_GW_INDUCED_SECOND_HELP_FUNC__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../general.h"
#include "../other.h"

//诱导引力波中与功率谱无关的积分函数

//方法二 1804.08577
int I_func_v_u_x(arb_t res, const arb_t v, const arb_t u, const arb_t x, slong prec);
int I_func_square_v_u_x_average(arb_t res, const arb_t v, const arb_t u, const arb_t x, slong prec);


#endif // __PBHS_GW_INDUCED_SECOND_HELP_FUNC__  


