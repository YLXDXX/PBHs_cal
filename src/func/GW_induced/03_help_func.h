#ifndef __PBHS_GW_INDUCED_THIRD_HELP_FUNC__   /* Include guard */
#define __PBHS_GW_INDUCED_THIRD_HELP_FUNC__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../general.h"
#include "../other.h"

//诱导引力波中与功率谱无关的积分函数

//方法三 2305.19950
void J_2_u_1_v_1_limited(arb_t res, const arb_t u_1, const arb_t v_1, slong prec);
void J_u_1_v_1_and_u_2_v_2_limited(arb_t res, const arb_t u_1, const arb_t v_1, const arb_t u_2, const arb_t v_2, slong prec);


#endif // __PBHS_GW_INDUCED_THIRD_HELP_FUNC__  


