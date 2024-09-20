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

// DOP853 方法，所需各系数
extern arb_t ODEs_DOP853_coe_c7, ODEs_DOP853_coe_c8, ODEs_DOP853_coe_c9, ODEs_DOP853_coe_c10, ODEs_DOP853_coe_c11, ODEs_DOP853_coe_c6, ODEs_DOP853_coe_c5, ODEs_DOP853_coe_c4, ODEs_DOP853_coe_c3, ODEs_DOP853_coe_c2, ODEs_DOP853_coe_b1, ODEs_DOP853_coe_b6, ODEs_DOP853_coe_b7, ODEs_DOP853_coe_b8, ODEs_DOP853_coe_b9, ODEs_DOP853_coe_b10, ODEs_DOP853_coe_b11;
extern arb_t ODEs_DOP853_coe_b12, ODEs_DOP853_coe_btilde1, ODEs_DOP853_coe_btilde6, ODEs_DOP853_coe_btilde7, ODEs_DOP853_coe_btilde8, ODEs_DOP853_coe_btilde9, ODEs_DOP853_coe_btilde10, ODEs_DOP853_coe_btilde11;
extern arb_t ODEs_DOP853_coe_btilde12, ODEs_DOP853_coe_er1, ODEs_DOP853_coe_er6, ODEs_DOP853_coe_er7, ODEs_DOP853_coe_er8, ODEs_DOP853_coe_er9, ODEs_DOP853_coe_er10, ODEs_DOP853_coe_er11, ODEs_DOP853_coe_er12, ODEs_DOP853_coe_a0201, ODEs_DOP853_coe_a0301;
extern arb_t ODEs_DOP853_coe_a0302, ODEs_DOP853_coe_a0401, ODEs_DOP853_coe_a0403, ODEs_DOP853_coe_a0501, ODEs_DOP853_coe_a0503, ODEs_DOP853_coe_a0504, ODEs_DOP853_coe_a0601, ODEs_DOP853_coe_a0604, ODEs_DOP853_coe_a0605, ODEs_DOP853_coe_a0701;
extern arb_t ODEs_DOP853_coe_a0704, ODEs_DOP853_coe_a0705, ODEs_DOP853_coe_a0706, ODEs_DOP853_coe_a0801, ODEs_DOP853_coe_a0804, ODEs_DOP853_coe_a0805, ODEs_DOP853_coe_a0806, ODEs_DOP853_coe_a0807, ODEs_DOP853_coe_a0901, ODEs_DOP853_coe_a0904;
extern arb_t ODEs_DOP853_coe_a0905, ODEs_DOP853_coe_a0906, ODEs_DOP853_coe_a0907, ODEs_DOP853_coe_a0908, ODEs_DOP853_coe_a1001, ODEs_DOP853_coe_a1004, ODEs_DOP853_coe_a1005, ODEs_DOP853_coe_a1006, ODEs_DOP853_coe_a1007, ODEs_DOP853_coe_a1008;
extern arb_t ODEs_DOP853_coe_a1009, ODEs_DOP853_coe_a1101, ODEs_DOP853_coe_a1104, ODEs_DOP853_coe_a1105, ODEs_DOP853_coe_a1106, ODEs_DOP853_coe_a1107, ODEs_DOP853_coe_a1108, ODEs_DOP853_coe_a1109, ODEs_DOP853_coe_a1110, ODEs_DOP853_coe_a1201;
extern arb_t ODEs_DOP853_coe_a1204, ODEs_DOP853_coe_a1205, ODEs_DOP853_coe_a1206, ODEs_DOP853_coe_a1207, ODEs_DOP853_coe_a1208, ODEs_DOP853_coe_a1209, ODEs_DOP853_coe_a1210, ODEs_DOP853_coe_a1211, ODEs_DOP853_coe_c14, ODEs_DOP853_coe_c15, ODEs_DOP853_coe_c16;
extern arb_t ODEs_DOP853_coe_a1401, ODEs_DOP853_coe_a1407, ODEs_DOP853_coe_a1408, ODEs_DOP853_coe_a1409, ODEs_DOP853_coe_a1410, ODEs_DOP853_coe_a1411, ODEs_DOP853_coe_a1412, ODEs_DOP853_coe_a1413, ODEs_DOP853_coe_a1501, ODEs_DOP853_coe_a1506;
extern arb_t ODEs_DOP853_coe_a1507, ODEs_DOP853_coe_a1508, ODEs_DOP853_coe_a1511, ODEs_DOP853_coe_a1512, ODEs_DOP853_coe_a1513, ODEs_DOP853_coe_a1514, ODEs_DOP853_coe_a1601, ODEs_DOP853_coe_a1606, ODEs_DOP853_coe_a1607, ODEs_DOP853_coe_a1608;
extern arb_t ODEs_DOP853_coe_a1609, ODEs_DOP853_coe_a1613, ODEs_DOP853_coe_a1614, ODEs_DOP853_coe_a1615, ODEs_DOP853_coe_d401, ODEs_DOP853_coe_d406, ODEs_DOP853_coe_d407, ODEs_DOP853_coe_d408, ODEs_DOP853_coe_d409, ODEs_DOP853_coe_d410, ODEs_DOP853_coe_d411;
extern arb_t ODEs_DOP853_coe_d412, ODEs_DOP853_coe_d413, ODEs_DOP853_coe_d414, ODEs_DOP853_coe_d415, ODEs_DOP853_coe_d416, ODEs_DOP853_coe_d501, ODEs_DOP853_coe_d506, ODEs_DOP853_coe_d507, ODEs_DOP853_coe_d508, ODEs_DOP853_coe_d509, ODEs_DOP853_coe_d510, ODEs_DOP853_coe_d511;
extern arb_t ODEs_DOP853_coe_d512, ODEs_DOP853_coe_d513, ODEs_DOP853_coe_d514, ODEs_DOP853_coe_d515, ODEs_DOP853_coe_d516, ODEs_DOP853_coe_d601, ODEs_DOP853_coe_d606, ODEs_DOP853_coe_d607, ODEs_DOP853_coe_d608, ODEs_DOP853_coe_d609, ODEs_DOP853_coe_d610, ODEs_DOP853_coe_d611;
extern arb_t ODEs_DOP853_coe_d612, ODEs_DOP853_coe_d613, ODEs_DOP853_coe_d614, ODEs_DOP853_coe_d615, ODEs_DOP853_coe_d616, ODEs_DOP853_coe_d701, ODEs_DOP853_coe_d706, ODEs_DOP853_coe_d707, ODEs_DOP853_coe_d708, ODEs_DOP853_coe_d709, ODEs_DOP853_coe_d710, ODEs_DOP853_coe_d711;
extern arb_t ODEs_DOP853_coe_d712, ODEs_DOP853_coe_d713, ODEs_DOP853_coe_d714, ODEs_DOP853_coe_d715, ODEs_DOP853_coe_d716;
extern arb_t  ODEs_DOP853_coe_bhh1,ODEs_DOP853_coe_bhh2,ODEs_DOP853_coe_bhh3;


//系数初始化
void ODEs_get_step_cal_coe(slong prec);
void ODEs_get_RFK45_cal_coe(slong prec);
void ODEs_get_DOPRI54_cal_coe(slong prec);
void ODEs_get_DOP853_cal_coe(slong prec);

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

ODEs_DOP853_dense_t ODEs_DOP853_dense_init(slong num, slong dim);
void ODEs_DOP853_dense_clear(ODEs_DOP853_dense_t dense_out, slong num, slong dim);

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

//常微分方程求解 DOP853 方法
void ODEs_DOP853(arb_ptr y_end, my_odes_func func, const slong dim, void *param, const slong order, //常微分方程组函数
                 const arb_t x_start, const arb_ptr y_start, //给定初始条件
                 const arb_t x_end, //求出点 x_end 对应的函数值
                 const arb_t error_abs, const arb_t error_rel, //绝对误差和相对误差
                 const ODEs_DOP853_dense_t dense_out, //此参数可以为空
                 slong prec);

#endif // __PBHS_MATH_METHOD_ODES__ 

