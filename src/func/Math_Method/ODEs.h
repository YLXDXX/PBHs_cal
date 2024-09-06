#ifndef __PBHS_MATH_METHOD_ODES__   /* Include guard */
#define __PBHS_MATH_METHOD_ODES__

#include <arb.h> //高精度实数运算
#include "../new_type.h"
#include "../phy_constant.h"
#include "func_constant.h"

//常微分方程求解

int ODEs_RFK45(arb_ptr res, my_odes_func func, const slong dim, void *param, const slong order, //常微分方程组函数
               const arb_t x_start, const arb_ptr y_start, //给定初始条件
               const arb_t x_end, //求出点 x_end 对应的函数值
               const slong num, const arb_t error, //num为迭代区间为[x_start,x_end]，将其分为几等份，从而给出初始步长
                                                   //误差为相对精度
               slong prec);

#endif // __PBHS_MATH_METHOD_ODES__ 

