#ifndef __PBHS_MATH_METHOD_FIT_DATE__   /* Include guard */
#define __PBHS_MATH_METHOD_FIT_DATE__ 

#include <arb.h> //高精度实数运算
#include "../new_type.h"
#include "../phy_constant.h"


//函数对 y=f(x) 将点 (x,y) 输出到文件，储存备用
int Func_output_point(int (*func)(arb_t f_res, const arb_t f_x, const slong f_order, slong prec),
                      const slong order,
                      const arb_t a, const arb_t b, const slong N, char* file, slong prec);

//数据拟合
int Func_fit_point(char* fitted_file, char* output_file,const unsigned int num_threads, slong prec);


//拟合数据拟合恢复
int Func_fit_restore(struct FUNC_FITTED_DATE **fit_res_2, char* fitted_file, slong prec);


//求相应拟合函数的值
int Fit_get_value(arb_t res, const arb_t x, const struct FUNC_FITTED_DATE* Fit_func,
                  const int order, slong prec);



#endif // __PBHS_MATH_METHOD_FIT_DATE__ 
