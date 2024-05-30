#ifndef __PBHS_GW_INDUCED_GW_FUNC__   /* Include guard */
#define __PBHS_GW_INDUCED_GW_FUNC__ 

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../general.h"
#include "../other.h"

//引力波产生时的能量密度与现今能量密度转换系数 c_g*Ω_{r,0}
void GW_spectra_convert_coefficient(arb_t res, slong prec);

#endif // __PBHS_GW_INDUCED_GW_FUNC__   
