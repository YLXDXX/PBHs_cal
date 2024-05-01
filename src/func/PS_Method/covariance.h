#ifndef __PBHS_PS_METHOD_COVARIANCE__   /* Include guard */
#define __PBHS_PS_METHOD_COVARIANCE__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../general.h"

int Variance_XX(arb_t res, const arb_t r, slong prec); //计算方差 XX
int Variance_XY(arb_t res, const arb_t r, slong prec); //计算方差 XY
int Variance_YY(arb_t res, const arb_t r, slong prec); //计算方差 YY


#endif // __PBHS_PS_METHOD_COVARIANCE__  
 
