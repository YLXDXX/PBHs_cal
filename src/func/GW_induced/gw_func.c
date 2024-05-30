#include "gw_func.h"

//引力波产生时的能量密度与现今能量密度转换系数 c_g*Ω_{r,0}
void GW_spectra_convert_coefficient(arb_t res, slong prec)
{
    arb_t s,t,w;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //转换系数 c_g * Ω_{r,0}
    arb_div(s,effective_g_star,effective_g_star_current,prec); //c_g= g_star/g_{star,0} * (g_{star,s,0}/g_{star,s})^{4/3}
    arb_div(t,effective_g_star_current_entropy,effective_g_star_entropy,prec);
    arb_one(w);
    arb_mul_ui(w,w,4,prec);
    arb_div_ui(w,w,3,prec);
    arb_pow(t,t,w,prec);
    arb_mul(w,s,t,prec);
    
    arb_mul(res,w,Omega_radiation,prec); //c_g * Ω_{r,0}
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
} 
