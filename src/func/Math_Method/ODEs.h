#ifndef __PBHS_MATH_METHOD_ODES__   /* Include guard */
#define __PBHS_MATH_METHOD_ODES__

#include <arb.h> //高精度实数运算
#include "../new_type.h"
#include "../phy_constant.h"
#include "func_constant.h"

//常微分方程求解

int ODEs_RFK45(arb_ptr y_end, my_odes_func func, const slong dim, void *param, const slong order, //常微分方程组函数
               const arb_t x_start, const arb_ptr y_start, //给定初始条件
               const arb_t x_end, //求出点 x_end 对应的函数值
               const slong num, const arb_t error, //num为迭代区间为[x_start,x_end]，将其分为几等份，从而给出初始步长
                                                   //误差为相对精度
               slong prec);



ODEs_point_output_t ODEs_point_output_init(slong num, slong dim); //初始化微分方程多点输出结构，为其分配内存
void ODEs_point_output_clear(ODEs_point_output_t p_out, slong num, slong dim); //清理微分方程多点输出结构，释放内存

//在给定区间 [x_start, x_end] 中，求解 n 个点，输出这 n 个点的坐标(x,y)，每个 x 对应多个 y
//这 n 个点包括区间端点，共把区间分为 n-1 份
int ODEs_RFK45_interval_point_output(ODEs_point_output_t p_out, 
                                     my_odes_func func, const slong dim, void *param, const slong order, //常微分方程组函数
                                     const arb_t x_start, const arb_ptr y_start, //给定初始条件
                                     const arb_t x_end, //给定区间 [x_start, x_end]
                                     const slong N, const arb_t error, //N为输出点个数，区间 N-1 等分，误差为相对精度
                                     slong prec);

#endif // __PBHS_MATH_METHOD_ODES__ 

