#ifndef __PBHS_MATH_METHOD_QUADRATURE_BINARY__   /* Include guard */
#define __PBHS_MATH_METHOD_QUADRATURE_BINARY__

#include <arb.h> //高精度实数运算
#include "../new_type.h"
#include "../phy_constant.h"
#include "func_constant.h"

#include "quadrature_binary_rectangle_gauss_kronrod.h"
#include "quadrature_binary_rectangle_double_exponential.h"

#include "quadrature_binary_func_gauss_kronrod.h"
#include "quadrature_binary_func_double_exponential.h"

//二元积分，积分区域为矩形
int integration_binary_rectangle(arb_t res, my_calc_func_binary func, void *param, const slong order,
                                 const arb_t x_a, const arb_t x_b, const arb_t x_error, 
                                 slong x_step_min , slong x_step_max,
                                 const arb_t y_a, const arb_t y_b, const arb_t y_error,
                                 slong y_step_min , slong y_step_max,
                                 slong prec);

//二元积分，积分区域非矩形
int integration_binary_func(arb_t res, my_calc_func_binary func, void* param, const slong order,
                            const arb_t x_a, const arb_t x_b, const arb_t x_error, 
                            slong x_step_min , slong x_step_max,
                            my_calc_func y_a_func, void *param_y_a,  const slong order_y_a,
                            my_calc_func y_b_func, void *param_y_b,  const slong order_y_b,
                            const arb_t y_error, slong y_step_min, slong y_step_max,
                            slong prec);

#endif // __PBHS_MATH_METHOD_QUADRATURE_BINARY__  
