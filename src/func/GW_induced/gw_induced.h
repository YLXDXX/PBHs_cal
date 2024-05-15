#ifndef __PBHS_GW_INDUCED_GW_INDUCED__   /* Include guard */
#define __PBHS_GW_INDUCED_GW_INDUCED__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../general.h"
#include "../other.h"

#include "01_power_spectra.h"
#include "02_power_spectra.h"
#include "03_power_spectra.h"


int GW_power_spectra(arb_t res, const arb_t eta, const arb_t k, slong prec); //诱导引力波的功率谱
int GW_current_energy_density(arb_t res, const arb_t k, slong prec); //当前的引力波能量密度



#endif // __PBHS_GW_INDUCED_GW_INDUCED__   
