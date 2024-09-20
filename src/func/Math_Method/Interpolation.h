#ifndef __PBHS_MATH_METHOD_INTERPOLATION__   /* Include guard */
#define __PBHS_MATH_METHOD_INTERPOLATION__ 

#include <arb.h> //高精度实数运算
#include "../new_type.h"
#include "../phy_constant.h"

//函数插值

//插值点位置查找
slong Interpolation_position_i_search(const arb_t x, const arb_ptr x_array, const slong N, slong prec );

Interp_coe_t Interpolation_coe_init(slong num); //初始化插值多项式系数，为其分配内存
void Interpolation_coe_clear(Interp_coe_t coe, slong num); //清理初始化插值多项式系数，释放内存

//通过各个点的离散函数值，给出相应的拟合函数 y=func(x)
void Interpolation_fit_func(arb_t res, const arb_t x,
                            const arb_ptr x_array, const arb_ptr y_array, const Interp_coe_t coe, const slong N,
                            slong prec);

//常微分方程求解中，RFK45 方法对应的内部插值法
void Interpolation_fit_func_odes_RFK45(arb_t res, const arb_t x,
                                       ODEs_RFK45_dense_t dense_out, const slong i_y, //i_y 表示微分方程解 y 中的第 i 个
                                       slong prec);


//常微分方程求解中，DOPRI54 方法对应的内部插值法
void Interpolation_fit_func_odes_DOPRI54(arb_t res, const arb_t x,
                                         ODEs_DOPRI54_dense_t dense_out, const slong i_y, //i_y 表示微分方程解 y 中的第 i 个
                                         slong prec);

//常微分方程求解中，DOP853 方法对应的内部插值法
void Interpolation_fit_func_odes_DOP853(arb_t res, const arb_t x,
                                        ODEs_DOP853_dense_t dense_out, const slong i_y, //i_y 表示微分方程解 y 中的第 i 个
                                        slong prec);

#endif // __PBHS_MATH_METHOD_INTERPOLATION__ 
