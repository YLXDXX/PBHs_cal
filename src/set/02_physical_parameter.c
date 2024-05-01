#include "header/02_physical_parameter.h"

void Set_physical_parameter(slong prec)
{
    //宇宙学基本常数
    //cold dark matter: Ω_DMC ≈ 0.27
    // baryons + dark matter: Ω_m ≈ 0.32
    //dark energy: Ω_Λ ≈ 0.68
    arb_set_str(Omega_DM,"0.27",prec); //暗物质所占比例
    arb_set_str(Omega_M,"0.32",prec); //baryons + dark matter 所占比例
    
    
    //哈勃常数 h
    //h = 0.740 ± 0.014 (supernovae)
    //h = 0.674 ± 0.005 (CMB)
    //H_0 ≡ 100*h kms^{−1} Mpc^{−1}
    arb_set_str(Hubble_constant_h,"0.674",prec);//采用CMB的观测
    
    //哈勃常数 H_0
    arb_mul_si(Hubble_constant_H_0, Hubble_constant_h, 100, prec);
    
    // 尺度因子 a(t)，某些公式的需要
    arb_one(Scale_factor_a);
    
    //尺度因子 a(t_0)
    arb_one(Scale_factor_a_0);
    
    //状态方程 w
    arb_one(Equation_Of_State_w); // w=1/3
    arb_div_ui(Equation_Of_State_w,Equation_Of_State_w,3,prec);
    
    
    arb_set_str(effective_g_star,"106.75",prec); //重新进入视界后，对应的有效自由度数目 g_*
    arb_set_str(effective_g_star_eq,"3.38",prec); //物质辐射相等的时刻，对应的有效自由度数目 g_{*,eq}
    arb_set_str(K_scale_eq,"0.07",prec); //在视界相等的时刻，对应的参考尺度 k_eq， K_scale_eq=0.07*Ω_M*h^2 
    arb_mul(K_scale_eq,K_scale_eq,Omega_M,prec);
    arb_mul(K_scale_eq,K_scale_eq,Hubble_constant_h,prec);
    arb_mul(K_scale_eq,K_scale_eq,Hubble_constant_h,prec);
    
}
