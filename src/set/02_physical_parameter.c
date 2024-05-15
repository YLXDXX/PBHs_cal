#include "header/02_physical_parameter.h"

void Set_physical_parameter(slong prec)
{
    //宇宙学基本常数
    
    //哈勃常数 h
    //h = 0.740 ± 0.014 (supernovae)
    //h = 0.674 ± 0.005 (CMB)
    //H_0 ≡ 100*h kms^{−1} Mpc^{−1}
    arb_set_str(Hubble_constant_h,"0.674",prec);//采用CMB的观测
    
    //哈勃常数 H_0
    arb_mul_si(Hubble_constant_H_0, Hubble_constant_h, 100, prec);
    
    
    //cold dark matter: Ω_DMC ≈ 0.27
    // baryons + dark matter: Ω_m ≈ 0.32
    //dark energy: Ω_Λ ≈ 0.68
    arb_set_str(Omega_DM,"0.27",prec); //暗物质所占比例
    arb_set_str(Omega_M,"0.32",prec); //baryons + dark matter 所占比例
    arb_set_str(Omega_radiation,"8.99E-5",prec); //这里采取的是Daniel Baumann's Cosmology Table 2.1 中的结果
                                                //辐射所占比例 https://pdg.lbl.gov/2019/reviews/rpp2019-rev-cosmological-parameters.pdf
                                                //https://arxiv.org/pdf/2305.19950
    arb_div(Omega_radiation,Omega_radiation,Hubble_constant_h,prec); //Ω_r=2.47*10^(−5)*h^(−2)
    arb_div(Omega_radiation,Omega_radiation,Hubble_constant_h,prec);
    
    
    // 尺度因子 a(t)，某些公式的需要
    arb_one(Scale_factor_a);
    
    //尺度因子 a(t_0)
    arb_one(Scale_factor_a_0);
    
    //状态方程 w
    arb_one(Equation_Of_State_w); // w=1/3
    arb_div_ui(Equation_Of_State_w,Equation_Of_State_w,3,prec);
    
    
    arb_set_str(effective_g_star,"106.75",prec); //重新进入视界后，对应的有效自由度数目 g_*
    arb_set_str(effective_g_star_eq,"3.38",prec); //物质辐射相等的时刻，对应的有效自由度数目 g_{*,eq}
    arb_set_str(effective_g_star_current,"3.363",prec); //参考 1609.04979
    arb_set_str(effective_g_star_current_entropy,"3.909",prec);
    
    arb_set_str(K_scale_eq,"0.07",prec); //在视界相等的时刻，对应的参考尺度 k_eq， K_scale_eq=0.07*Ω_M*h^2 
    arb_mul(K_scale_eq,K_scale_eq,Omega_M,prec);
    arb_mul(K_scale_eq,K_scale_eq,Hubble_constant_h,prec);
    arb_mul(K_scale_eq,K_scale_eq,Hubble_constant_h,prec);
    
}
