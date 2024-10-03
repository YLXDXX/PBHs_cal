#ifndef __PBHS_GENERAL_NON_GUASSIANITY__   /* Include guard */
#define __PBHS_GENERAL_NON_GUASSIANITY__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"

#include "basis.h"


int Non_Gaussianity_exponential_tail_n(arb_t res, const arb_t x, const slong order, slong prec); // exponential_tail 类型函数及其导数

int Non_Gaussianity_up_step_n(arb_t res, const arb_t x, const slong order, slong prec); // up_step 类型函数及其导数

int Non_Gaussianity_power_expansion_n(arb_t res, const arb_t x, const slong order, slong prec); // power_expansion 类型函数及其导数

int Non_Gaussianity_narrow_1_2_up_step_n(arb_t res, const arb_t x, const slong order, slong prec); //有限宽 up_step, 扰动 1+2 阶
int Non_Gaussianity_narrow_1_up_step_n(arb_t res, const arb_t x, const slong order, slong prec); //有限宽 up_step, 扰动 1 阶


#endif // __PBHS_GENERAL_NON_GUASSIANITY__ 

