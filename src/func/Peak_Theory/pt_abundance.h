#ifndef __PBHS_PEAK_THEORY_ABUNDANCE__   /* Include guard */
#define __PBHS_PEAK_THEORY_ABUNDANCE__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../general.h"


#include "pt_generate_mass.h"
#include "number_density.h"

int PT_abundance_beta_m(arb_t res, const arb_t M, slong prec); //PBH abundance beta_PBH(M)
int PT_abundance_beta_all(arb_t res, slong prec); //PBH abundance beta_all

int PT_abundance_f_m(arb_t res, const arb_t M, slong prec); //PBH abundance f_PBH(M)
int PT_abundance_f_all(arb_t res, slong prec); //PBH abundance f_all


#endif // __PBHS_PEAK_THEORY_ABUNDANCE__
