#ifndef __PBHS_MATH_METHOD_QUADRATURE_BINARY_RECTANGLE_GAUSS_KRONROD__   /* Include guard */
#define __PBHS_MATH_METHOD_QUADRATURE_BINARY_RECTANGLE_GAUSS_KRONROD__

#include <arb.h> //高精度实数运算
#include "../new_type.h"
#include "../phy_constant.h"
#include "func_constant.h"

//二元积分，积分区域为矩形
int integration_binary_rectangle_gauss_kronrod(arb_t res, my_calc_func_binary func, void *param, const slong order,
                                               const arb_t x_a, const arb_t x_b, const arb_t x_error, 
                                               slong x_step_min , slong x_step_max,
                                               const arb_t y_a, const arb_t y_b, const arb_t y_error,
                                               slong y_step_min , slong y_step_max,
                                               slong prec);


#endif // __PBHS_MATH_METHOD_QUADRATURE_BINARY_RECTANGLE_GAUSS_KRONROD__  
