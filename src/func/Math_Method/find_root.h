#ifndef __PBHS_MATH_METHOD_FIND_ROOT__   /* Include guard */
#define __PBHS_MATH_METHOD_FIND_ROOT__

#include <stdbool.h> // bool 类型
#include <arb.h> //高精度实数运算
#include "../new_type.h"

extern bool Stdout_verbose; //命令行输出

extern enum FIND_ROOT_ENUM_TYPE Find_root_method; //找根方法


//Brent’s method 求根
int Find_root_brent_method(arb_t res, my_calc_func func, void *param, const slong order,
                           const arb_t a, const arb_t b, const arb_t error,
                           slong prec);

//在指定的区间内找根算法
//指定函数，指定区间，指定子区间数，指定精度
int Find_interval_root(arb_t res, my_calc_func func, void *param, const slong order,
                       const arb_t a, const arb_t b, const arb_t error,
                       const slong num, const enum Root_type r_type, slong prec);


//在一个给定的范围内，尽可能的找出该范围内的所有根
int Find_interval_multi_root(arb_ptr* res, my_calc_func func, void *param, const slong order,
                             const arb_t a, const arb_t b, const arb_t error,
                             const slong num, slong prec);


#endif // __PBHS_MATH_METHOD_FIND_ROOT__ 
