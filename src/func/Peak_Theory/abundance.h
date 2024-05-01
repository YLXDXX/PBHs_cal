#ifndef __PBHS_PEAK_THEORY_ABUNDANCE__   /* Include guard */
#define __PBHS_PEAK_THEORY_ABUNDANCE__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../general.h"


#include "generate_mass.h"
#include "number_density.h"

int PBH_abundance_f_to_M(arb_t res,const arb_t M, slong prec); //PBH abundance f_PBH(M)

#endif // __PBHS_PEAK_THEORY_ABUNDANCE__
