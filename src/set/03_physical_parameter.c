#include "header/03_physical_parameter.h"

void Set_physical_parameter(slong prec)
{
    //宇宙学基本常数
    
    //若 Planck 2018 有相关的值，就采用 Planck 2018 的值 1807.06209
    //Page 16, TT,TE,EE+lowE+lensing+BAO 68% limits
    
    //哈勃常数 h
    //h = 0.740 ± 0.014 (supernovae)
    //h = 0.674 ± 0.005 (CMB)
    //H_0 ≡ 100*h kms^{−1} Mpc^{−1}
    arb_set_str(Hubble_constant_h,"0.6766",prec);//采用CMB的观测
    
    //哈勃常数 H_0
    arb_mul_si(Hubble_constant_H_0, Hubble_constant_h, 100, prec);
    
    
    //cold dark matter: Ω_DMC ≈ 0.27
    // baryons + dark matter: Ω_m ≈ 0.32
    //dark energy: Ω_Λ ≈ 0.68
    arb_set_str(Omega_DM,"0.11933",prec); //暗物质所占比例 cold dark matter density //Ω_c*h^2=0.11933
    arb_div(Omega_DM,Omega_DM,Hubble_constant_h,prec);
    arb_div(Omega_DM,Omega_DM,Hubble_constant_h,prec);
    
    arb_set_str(Omega_B,"0.02242",prec); //普通重子物质所占比例 baryons //Ω_b*h^2=0.02242
    arb_div(Omega_B,Omega_B,Hubble_constant_h,prec);
    arb_div(Omega_B,Omega_B,Hubble_constant_h,prec);
    
    arb_set_str(Omega_M,"0.3111",prec); //baryons + dark matter 所占比例
    
    arb_set_str(Omega_radiation,"8.99E-5",prec); //辐射占比，radiation density
                                                //这里采取的是Daniel Baumann's Cosmology Table 2.1 中的结果
                                                //辐射所占比例 https://pdg.lbl.gov/2019/reviews/rpp2019-rev-cosmological-parameters.pdf
                                                //Ω_γ=2.47*10^{−5}*h^(−2)
                                                //https://arxiv.org/pdf/2305.19950 --> Ω_{rad,0} = 4.2*10^{-5}*h^(−2)
    //辐射除了光子外，还有中微子
    arb_set_str(Omega_gamma,"5.35E-5",prec); //光子所占的比例 photon density
    
    arb_set_str(Omega_Lambda,"0.6889",prec); //暗能量 Ω_Λ 所占比例
    
    // 尺度因子 a(t)，某些公式的需要
    arb_one(Scale_factor_a);
    
    //尺度因子 a(t_0)
    arb_one(Scale_factor_a_0);
    
    //状态方程 w
    arb_one(Equation_Of_State_w); // w=1/3
    arb_div_ui(Equation_Of_State_w,Equation_Of_State_w,3,prec);
    //注意到，在本程序的计算中，只有少数地方考虑了状态方程的影响
    //大部分地方，都没采用的 w=1/3 的值来表达计算的，故若要考虑 w≠1/3 的情况
    //并不能仅仅修改这里的值，还需要做大量函数的匹配
    //ToDO：对所有与状态方程相关的函数做变量 w 的匹配
    
    //物质辐射相等时期的一些量
    /*
    arb_set_str(K_scale_eq,"0.07",prec); //在视界相等的时刻，对应的参考尺度 k_eq， K_scale_eq=0.07*Ω_M*h^2
                                        //参考 1904.10298
    arb_mul(K_scale_eq,K_scale_eq,Omega_M,prec);
    arb_mul(K_scale_eq,K_scale_eq,Hubble_constant_h,prec);
    arb_mul(K_scale_eq,K_scale_eq,Hubble_constant_h,prec);
    */
    //参考 1812.00674, T_eq = (1 + z_eq)*2.725 K, z_eq ≃2.4*10^4*Ω_m*h^2
    //1+z=1/a --> T_eq = 2.725 k / a_eq
    //Planck 2018 给出了这些值的限制，直接采用 1807.06209
    
    arb_set_str(K_scale_eq,"0.010339",prec); //物质辐射相等时的共动波数，单位 Mpc^-1
    
    arb_set_str(Z_scale_eq,"3387",prec); //物质辐射相等时的红移
    
    arb_set(A_scale_eq,Z_scale_eq); //物质辐射相等时的尺度因子 a=1/(1+z)
    arb_add_ui(A_scale_eq,A_scale_eq,1,prec);
    arb_inv(A_scale_eq,A_scale_eq,prec);
    
    arb_set_str(T_scale_eq,"2.725",prec); //物质辐射相等时的温度 T_eq=(1+z_eq)*2.725 K = 2.725 K / a_eq, 1812.00674
    arb_div(T_scale_eq,T_scale_eq,A_scale_eq,prec);
    //这里温度的单位是 k 需转为 Gev, 1 k = 8.61732814974056E-05 eV
    arb_t T_scale_eq_temp;
    arb_init(T_scale_eq_temp);
    arb_set_str(T_scale_eq_temp,"8.61732814974056E-05",prec);
    arb_mul(T_scale_eq,T_scale_eq,T_scale_eq_temp,prec);
    arb_div_ui(T_scale_eq,T_scale_eq,1E9,prec); //1 Gev = 10^3 Mev = 10^9 ev
    arb_clear(T_scale_eq_temp);
    
    
    //相对论自由度数和熵自由度数设置
    //通过相应的拟合公式给出的 1803.01038
    arb_set_str(effective_g_star_eq,"3.383",prec); //物质辐射相等的时刻，对应的相对论有效自由度数目 g_{*,eq}
    arb_set_str(effective_g_star_eq_entropy,"3.931",prec); //物质辐射相等的时刻，对应的熵有效自由度数目 g_{*,s,eq}
    
    //当今的相对论自由度数和熵自由度数，等于物质辐射相等时刻的值
    arb_set_str(effective_g_star_current,"3.363",prec); //参考 1609.04979
    arb_set_str(effective_g_star_current_entropy,"3.909",prec);
    
    //由自由度数拟合函数，自动通过 k_star 得到 其对应的自由度数，当功率谱中心位移移动后，不用需手动调整
    //当前仅对于 delta, log-normal, broken power law 功率谱 有效，对于其它类型功率谱，需手动设置
    if( Power_spectrum_type==lognormal_type || Power_spectrum_type==delta_type || Power_spectrum_type==broken_power_law_type )
    {
        Func_k_to_degrees_of_freedom(effective_g_star, effective_g_star_entropy, K_star, prec);
    }else if ( Power_spectrum_type==upward_step_spectra_type || Power_spectrum_type==numerical_cal_type )
    {
        Func_k_to_degrees_of_freedom(effective_g_star, effective_g_star_entropy, Upward_step_spectra_k_c, prec);
    }
    else
    {
        printf("注意，检测到功率谱类型不为 delta | log-normal | broken power law | upward step, 对应的自由度数目需手动设置\n ");
        arb_set_str(effective_g_star,"106.75",prec); //重新进入视界后，对应的相对论有效自由度数目 g_*
        arb_set_str(effective_g_star_entropy,"106.75",prec); //重新进入视界后，对应的熵有效自由度数目 g_{*,s}
    }
    
}
