#ifndef __PBHS_OTHER_OTHER_O__   /* Include guard */
#define __PBHS_OTHER_OTHER_O__

#include <arb.h> //高精度实数运算
#include "../new_type.h"
#include "../phy_constant.h"

int Heaviside_Theta_function(arb_t res,const arb_t x,slong prec); //阶梯函数
void Func_GW_f_to_k(arb_t k, const arb_t f, slong prec); //诱引引力波频率 f nHZ 到波数 k Mpc^-1

//线性分割，类似于 numpy.linspace
void Get_interval_linspace_point(arb_ptr x, const arb_t a, const arb_t b, const slong N, slong prec);

//对数分割，类似于 numpy.logspace
void Get_interval_logspace_point(arb_ptr x, const arb_t a, const arb_t b, const slong N, slong prec);

#endif // __PBHS_OTHER_OTHER_O__ 
