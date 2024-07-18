#ifndef __PBHS_GENERAL_TYPICAL_PROFILE__   /* Include guard */
#define __PBHS_GENERAL_TYPICAL_PROFILE__

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"

#include "basis.h"
#include "non_gaussianity.h"

int Help_sigma_n_square(arb_t res, const slong n,slong prec); //计算 gradient moments σ_n

int Help_sinc_n(arb_t res, const arb_t x, const slong order, slong prec); //计算 sin(x)/x 及其各阶导数

int Help_psi_n(arb_t res, const arb_t r, const slong order, slong prec); //计算 ψ_1(r)

int zeta_Gauss_profile_n(arb_t res, const arb_t r, const slong order, slong prec); //ζ_G 及其各阶导数

int zeta_profile_n(arb_t res, const arb_t r, const slong order, slong prec); // ζ(r) 及其各阶导数


#endif // __PBHS_GENERAL_TYPICAL_PROFILE__ 
