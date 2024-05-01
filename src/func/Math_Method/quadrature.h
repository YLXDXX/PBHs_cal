#ifndef __PBHS_MATH_METHOD_QUADRATURE__   /* Include guard */
#define __PBHS_MATH_METHOD_QUADRATURE__

#include <arb.h> //高精度实数运算
#include "../new_type.h"
#include "../phy_constant.h"
#include "func_constant.h"

//积分中的排序算法，暴露出来给 quadrature_binaty.c 中的二元函数积分用
void small2big_order(arb_ptr as, arb_ptr bs, arb_ptr es, arb_ptr gs, slong n);


//一元函数积分
//Simpson公式
int simpson_formula(arb_t res, int (*func)(arb_t f_res, const arb_t x, slong prec),
                    const arb_t a, const arb_t b, slong prec); //Simpson公式

//Simpsont积分
int integration_simpson(arb_t int_sum, int (*func)(arb_t f_res, const arb_t x, slong prec),
                        const arb_t a, const arb_t b,const arb_t ans,const arb_t eps,
                        slong step_min , slong step_max,
                        slong prec);  //Simpsont积分，加入最少迭代次数，加入最大迭代次数

//Gauss–Kronrod积分法，自适应，15/65点的
//Gauss–Kronrodt积分适用于通常的和有振荡行为的函数
//Gauss–Kronrodt积分不适用于区间端点有发散的情况，及无穷区间的情况
int integration_gauss_kronrod_iterate(arb_t res, my_calc_func func, void *param, const slong order,
                              const arb_t a, const arb_t b, const arb_t error,
                              slong step_min , slong step_max,
                              slong prec); //自适应计算

int integration_gauss_kronrod_recursive(arb_t res, my_calc_func func, void *param, const slong order,
                                      const arb_t a, const arb_t b, const arb_t error,
                                      slong step_min , slong step_max,
                                      slong prec); //自适应计算

//Double Exponential积分，可适用于无穷区间积分，及区间端点处函数发散的情况
//特别是在区间端点处函数发散的情况，Double Exponential积分明显好于Gauss–Kronrod积分
int Double_Exponential_Quadrature(arb_t res, my_calc_func func, void *param, const slong order,
                                  const arb_t a, const arb_t b, const arb_t eps,
                                  slong step_min , slong step_max,
                                  slong prec);

int Integration_arb(arb_t res, my_calc_func func, void *param, const slong order,
                    const arb_t a, const arb_t b, const arb_t error,
                    slong step_min , slong step_max,
                    slong prec);

#endif // __PBHS_MATH_METHOD_QUADRATURE__ 

