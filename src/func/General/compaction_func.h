#ifndef __PBHS_GENERAL_COMPACTION_FUNC__   /* Include guard */
#define __PBHS_GENERAL_COMPACTION_FUNC__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"

#include "basis.h"
#include "non_gaussianity.h"
#include "typical_profile.h"

void interior_assist_f_w(arb_t res, const arb_t w, slong prec); //f(w)

int C_r_profile_n(arb_t res, const arb_t r, const slong order, slong prec); // C(r) 及其各阶导数

int Find_r_max(arb_t res, slong prec); //求C(r)的最大值，在区间[a, b]中找


#endif // __PBHS_GENERAL_COMPACTION_FUNC__ 
