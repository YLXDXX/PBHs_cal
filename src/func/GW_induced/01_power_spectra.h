#ifndef __PBHS_GW_INDUCED_FIRST_POWER_SPECTRA__   /* Include guard */
#define __PBHS_GW_INDUCED_FIRST_POWER_SPECTRA__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../general.h"
#include "../other.h"

#include "01_help_func.h"

//方法一 1804.07732
int GW_power_spectra_01(arb_t res, const arb_t eta, const arb_t k, slong prec); //诱导引力波的功率谱
int GW_current_energy_density_01(arb_t res, const arb_t k, slong prec); //当前的引力波能量密度



#endif // __PBHS_GW_INDUCED_FIRST_POWER_SPECTRA__  


