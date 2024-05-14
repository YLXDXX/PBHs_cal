#ifndef __PBHS_MATH_METHOD_QUADRATURE_BINARY_ADAPTIVE__   /* Include guard */
#define __PBHS_MATH_METHOD_QUADRATURE_BINARY_ADAPTIVE__

#include <arb.h> //高精度实数运算
#include "../new_type.h"
#include "../phy_constant.h"
#include "func_constant.h"

#include "quadrature_binary_rectangle_gauss_kronrod.h"
#include "quadrature_binary_rectangle_double_exponential.h"

#include "quadrature_binary_func_gauss_kronrod.h"
#include "quadrature_binary_func_double_exponential.h"


void small2big_order_2D(arb_ptr x_as, arb_ptr x_bs, arb_ptr y_as, arb_ptr y_bs, arb_ptr x_es, arb_ptr y_es, arb_ptr gs, slong n);

//二元积分，对矩形积分区域采用自适应算法
int integration_binary_rectangle_adaptive_gauss_kronrod(arb_t res, my_calc_func_binary func, void *param, const slong order,
                                 const arb_t x_a, const arb_t x_b, const arb_t x_error, 
                                 slong x_step_min , slong x_step_max,
                                 const arb_t y_a, const arb_t y_b, const arb_t y_error,
                                 slong y_step_min , slong y_step_max,
                                 slong prec);

int integration_binary_rectangle_adaptive_double_exponential(arb_t res, my_calc_func_binary func, void *param, const slong order,
                                                             const arb_t x_a, const arb_t x_b, const arb_t x_error,
                                                             slong x_step_min , slong x_step_max,
                                                             const arb_t y_a, const arb_t y_b, const arb_t y_error,
                                                             slong y_step_min , slong y_step_max, 
                                                             slong prec);

#endif // __PBHS_MATH_METHOD_QUADRATURE_BINARY_ADAPTIVE__  
