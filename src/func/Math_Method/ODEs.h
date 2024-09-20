#ifndef __PBHS_MATH_METHOD_ODES__   /* Include guard */
#define __PBHS_MATH_METHOD_ODES__

#include <arb.h> //高精度实数运算
#include "../new_type.h"
#include "../phy_constant.h"
#include "func_constant.h"

//系数：估算误差，用于自动步长
extern arb_t ODEs_step_max_factor,ODEs_step_min_factor,ODEs_step_factor;

//方法 RFK45 系数
extern arb_ptr ODEs_RFK45_coe_ah;
extern arb_ptr ODEs_RFK45_coe_b1,ODEs_RFK45_coe_b2;
extern arb_ptr ODEs_RFK45_coe_b3,ODEs_RFK45_coe_b4,ODEs_RFK45_coe_b5;
extern arb_ptr ODEs_RFK45_coe_ck_s,ODEs_RFK45_coe_ck_t;

// DOPRI54 方法，所需各系数
extern arb_ptr ODEs_DOPRI54_coe_c;
extern arb_ptr ODEs_DOPRI54_coe_b_up;
extern arb_ptr ODEs_DOPRI54_coe_b_dw;
extern arb_ptr ODEs_DOPRI54_coe_a_1;
extern arb_ptr ODEs_DOPRI54_coe_a_2;
extern arb_ptr ODEs_DOPRI54_coe_a_3;
extern arb_ptr ODEs_DOPRI54_coe_a_4;
extern arb_ptr ODEs_DOPRI54_coe_a_5;
extern arb_ptr ODEs_DOPRI54_coe_a_6;
extern arb_ptr ODEs_DOPRI54_coe_dense_out;


//系数初始化
void ODEs_get_step_cal_coe(slong prec);
void ODEs_get_RFK45_cal_coe(slong prec);
void ODEs_get_DOPRI54_cal_coe(slong prec);

//计算平均误差，用于步长估算
void ODEs_RMS_norm(arb_t res, const arb_ptr x, const slong dim, const arb_ptr sc_i, slong prec);

//计算初始迭代步长
void ODEs_select_initial_step(arb_t h, my_odes_func func, const slong dim, void *param, const slong order, //常微分方程组函数
                              const arb_t x_start, const arb_ptr y_start, //给定初始条件
                              const arb_ptr yp_start, //func(x_start,y_start)
                              const arb_t gap_size, //计算区间间隔大小
                              const slong q, //估算方法的阶
                              int direction,
                              const arb_t error_abs, const arb_t error_rel, //绝对误差和相对误差 //误差为绝对精度
                              slong prec);

//为处理相关问题，定义了一些结构体，其相关操作
ODEs_RFK45_dense_t ODEs_RFK45_dense_init(slong num, slong dim);
void ODEs_RFK45_dense_clear(ODEs_RFK45_dense_t dense_out, slong num, slong dim);

ODEs_DOPRI54_dense_t ODEs_DOPRI54_dense_init(slong num, slong dim);
void ODEs_DOPRI54_dense_clear(ODEs_DOPRI54_dense_t dense_out, slong num, slong dim);


//常微分方程求解 RFK45 方法
void ODEs_RFK45(arb_ptr y_end, my_odes_func func, const slong dim, void *param, const slong order, //常微分方程组函数
               const arb_t x_start, const arb_ptr y_start, //给定初始条件
               const arb_t x_end, //求出点 x_end 对应的函数值
               const arb_t error_abs, const arb_t error_rel, //绝对误差和相对误差
               const ODEs_RFK45_dense_t dense_out, //此参数可以为空
               slong prec);

//常微分方程求解 DOPRI54 方法
void ODEs_DOPRI54(arb_ptr y_end, my_odes_func func, const slong dim, void *param, const slong order, //常微分方程组函数
                  const arb_t x_start, const arb_ptr y_start, //给定初始条件
                  const arb_t x_end, //求出点 x_end 对应的函数值
                  const arb_t error_abs, const arb_t error_rel, //绝对误差和相对误差 //误差为绝对精度
                  const ODEs_DOPRI54_dense_t dense_out, //此参数可以为空
                  slong prec);


#endif // __PBHS_MATH_METHOD_ODES__ 

