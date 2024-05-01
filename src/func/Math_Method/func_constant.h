#ifndef __PBHS_MATH_METHOD_FUNC__CONSTANT__   /* Include guard */
#define __PBHS_MATH_METHOD_FUNC__CONSTANT__

#include <arb.h> //高精度实数运算
#include "../new_type.h"

extern enum INTEGRAL_ENUM_TYPE Integral_method; //积分方法

extern arb_ptr INT_GUSS_KRONROD_COFFI; // gauss_kronrod 积分系数

extern arb_t INT_MIN_INTERVAL; //积分时，子区间最小间隔

int get_gauss_kronrod_node_weight(const int num, slong prec); //获取 gauss_kronrod 的节点位置和权重



#endif // __PBHS_MATH_METHOD_FUNC__CONSTANT__
