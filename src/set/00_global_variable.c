#include "header/00_global_variable.h" 

void Set_global_variable(slong prec)
{
    //积分时，子区间最小间隔
    arb_init(INT_MIN_INTERVAL); 
    
    //数学常数
    arb_init(Pi); //常数π 
    arb_init(Pi_2); //常数2π (πx2)
    
    //宇宙学基本常数
    arb_init(Omega_DM); //暗物质所占比例
    arb_init(Omega_M); //baryons + dark matter 所占比例
    arb_init(Omega_B); //普通重子物质 baryons
    arb_init(Omega_Lambda); //暗能量 Ω_Λ
    arb_init(Omega_radiation); //辐射所占的比例 radiation density
    arb_init(Omega_gamma); //光子所占的比例 photon density
    arb_init(Hubble_constant_h); //哈勃常数 h
    arb_init(Hubble_constant_H_0); //哈勃常数 H_0
    arb_init(Scale_factor_a);// 尺度因子 a(t)，某些公式的需要
    arb_init(Scale_factor_a_0); //尺度因子 a(t_0)
    arb_init(Equation_Of_State_w); //状态方程 w
    arb_init(K_scale_eq); //在物质辐射相等时刻，对应的共动波数 k_eq
    arb_init(Z_scale_eq); //在物质辐射相等时刻，对应的红移 z_eq
    arb_init(A_scale_eq); //在物质辐射相等时刻，对应的尺度因子 a_eq
    arb_init(T_scale_eq); //在物质辐射相等时刻，对应的温度 T_eq
    
    arb_init(effective_g_star_eq); //在物质辐射相等时刻，对应的相对论有效自由度数目
    arb_init(effective_g_star_eq_entropy); //在物质辐射相等时刻，对应的熵有效自由度数目
    arb_init(effective_g_star); //重新进入视界后，对应的相对论有效自由度数目
    arb_init(effective_g_star_entropy); //重新进入视界后，对应的熵有效自由度数目
    arb_init(effective_g_star_current); //当今对应的相对论有效自由度数目
    arb_init(effective_g_star_current_entropy); //当今对应的熵有效自由度数目
    
    
    //各功率谱的设定
    arb_init(Power_A); //功率谱振幅
    arb_init(Power_sigma); //功率谱掌宽，例如 log-normal 的
    arb_init(Log_normal_mul_sigma); //功率谱掌宽的倍数
    arb_init(K_star); //lognormal 功率谱参考尺度，单位 Mpc^{-1}
    arb_init(Ln_K_star); //功率谱参考尺度，单位 Mpc^{-1}, 取对数
    arb_init(Box_K_middle);//Box功率谱 中心点
    arb_init(Box_width); //Box功率谱 宽度
    arb_init(BPL_alpha); //broken power law功率谱使用 α
    arb_init(BPL_beta); //broken power law功率谱使用 β
    arb_init(BPL_gamma); //broken power law功率谱使用 γ
    arb_init(Link_CMB_K_t); //link cmb 功率谱用
    arb_init(Link_CMB_K_m);
    arb_init(Link_CMB_A_s);
    arb_init(Link_CMB_A_b);
    arb_init(Link_CMB_K_star);
    arb_init(Link_CMB_P_m);
    arb_init(Link_CMB_n_s);
    arb_init(Link_CMB_n_b);
    arb_init(Upward_step_spectra_epsilon_S); //upward step模型的原生功率谱, SR-USR-SR
    arb_init(Upward_step_spectra_epsilon_V);
    arb_init(Upward_step_spectra_eta_V);
    arb_init(Upward_step_spectra_g);
    arb_init(Upward_step_spectra_h);
    arb_init( Upward_step_spectra_V_0);
    arb_init(Upward_step_spectra_Delta_V);
    arb_init(Upward_step_spectra_Delta_phi_USR);
    arb_init(Upward_step_spectra_phi_s);
    arb_init(Upward_step_spectra_phi_c);
    arb_init(Upward_step_spectra_k_s);
    arb_init(Upward_step_spectra_k_c);
    arb_init(Upward_step_spectra_ln_k_s);
    arb_init(Upward_step_spectra_ln_k_c);
    arb_init(Upward_step_spectra_tau_c);
    
    //为求ζ(r)所需的积分区间设定
    arb_init(Int_sigma_n_min); //sigma_n ，初始化各变量
    arb_init(Int_sigma_n_max);
    arb_init(Int_sigma_n_precision);
    
    //计算ζ(r)的辅助变量
    arb_init(Sigma_0_square); //(σ_n)^2
    arb_init(Sigma_1_square);
    arb_init(Sigma_2_square);
    arb_init(Sigma_3_square);
    arb_init(Sigma_4_square);
    arb_init(Gamma_1); //(γ_1)
    arb_init(Gamma_3); //(γ_3)
    arb_init(R_1); //(R_1)
    arb_init(R_3); //(R_3)
    
    
    arb_init(R_MAX); // ζ(r) 取最大值时的 r 值
    arb_init(PT_mu); // ζ(r) 的参数 µ
    arb_init(PT_k); // ζ(r) 的参数 k
    arb_init(PT_beta_cal_need_exp_2_zeta_m); // 某个 M 对应的 e^(2*ζ_m)
    
    
    //主要计算参数
    arb_init(Int_r_min); // 求 r_m 的求值区间和精度
    arb_init(Int_r_max);
    arb_init(Int_r_precision); 
    
    arb_init(Int_mu_min); // 求 μ_th 的区间和精度
    arb_init(Int_mu_max);
    arb_init(Int_mu_precision); 
    
    arb_init(C_m_average_precision); //求 C_m_average 精度
    
    arb_init(Root_M_to_mu_min); //peak theory： M -> μ 求根用
    arb_init(Root_M_to_mu_max);
    arb_init(Root_M_to_mu_precision);
    
    arb_init(Int_n_pk_k_min); //peak theory： n_pk(mu,k) 中 k 的积分区间
    arb_init(Int_n_pk_k_max);
    arb_init(Int_n_pk_k_precision);
    
    arb_init(PS_Int_variance_min); // PS 计算方差 XX YY XY 积分区精和精度
    arb_init(PS_Int_variance_max);
    arb_init(PS_Int_variance_precision);
    
    arb_init(PS_Int_P_C_l_min); // PS 计算 C_ℓ 的概率密度分布 P(C_l) 积分区间和精度
    arb_init(PS_Int_P_C_l_max);
    arb_init(PS_Int_P_C_l_precision);
    
    arb_init(PS_Root_C_l_to_Y_min); // PS中δ情况，计算P(C_ℓ)需反解Y，求解区间和精度
    arb_init(PS_Root_C_l_to_Y_max);
    arb_init(PS_Root_C_l_to_Y_precision);
    
    arb_init(PS_abundance_f_all_precision); //PS 最终占比f积分的精度
    
    
    //非高斯相关参数设定
    //exponential_tail_type
    arb_init(Exponential_tail_beta);//exponential_tail_type 中β的取值，文献中通常取为β=3
    
    //power_expansion_type
    arb_init(Power_expansion_f); //power-series expansion 二次项 f_NL -> A
    arb_init(Power_expansion_g); //power-series expansion 三次项 g_NL -> B
    arb_init(Power_expansion_four); //power-series expansion 四次项 four -> C
    arb_init(Power_expansion_five); //power-series expansion 五次项 five -> D
    arb_init(Power_expansion_six); //power-series expansion 六次项 six -> E
    
    //up_step_type
    arb_init(Up_step_h); //up_step_type 其中Up_step_h取值为负
    
    
    arb_init(R_K_to_r_m); //求r_m动态区间用
    
    arb_init(PS_abundance_beta_delta_k_all_min); //用δ谱求连续谱，计算∫β(k_◦)dk_◦用，
    arb_init(PS_abundance_beta_delta_k_all_max);
    arb_init(PS_abundance_beta_delta_k_all_precision);
    
    arb_init(PS_abundance_beta_delta_k_M_precision); //考虑所有k计算β(M)用
    
    
    arb_init(PT_mu_th); // ζ(r) 的参数 µ 的临界值
    arb_init(PT_mu_max); // ζ(r) 的参数 µ 的临界值
    arb_init(Q_parameter_th); //ζ(r)取临界值时的q参数
    
    arb_init(Mass_K); //临界坍缩相关参数
    arb_init(Mass_gamma);
    
    arb_init(R_m_times_K); //对于δ谱，r_m*k_*为一定值
    
    
    arb_init(PS_Sigma_XX);//PS 计算方差 XX YY XY 的值
    arb_init(PS_Sigma_XY);
    arb_init(PS_Sigma_YX);
    arb_init(PS_Sigma_YY);
    arb_init(PS_Sigma_XX_save);
    arb_init(PS_Sigma_XY_save);
    arb_init(PS_Sigma_YY_save);
    
    arb_init(PS_C_th); //compaction function threshold
    arb_init(PS_C_l_th); // linear compaction function threshold
    arb_init(PS_M_ratio_max); //相对质量的最大比
    arb_init(PT_M_ratio_max);
    arb_init(PS_M_ratio_min); //相对质量的最小比
    arb_init(PT_M_ratio_min);
    
    arb_init(Continuum_spectrum_x_m); //利用δ谱的x_m，求连续谱的特征模式 k_ch
    arb_init(Continuum_spectrum_k_ch); //连续谱的特征模式
    arb_init(Continuum_spectrum_A_old);
    
    arb_init(Func_output_x_min);//函数数据的输入输出，及拟合
    arb_init(Func_output_x_max);
    
    //诱导引力波相关参数设定
    arb_init(Int_GW_I_func_min); // 变量初始化
    arb_init(Int_GW_I_func_max);
    arb_init(Int_GW_I_func_precision);
    arb_init(Int_GW_power_spectra_min);
    arb_init(Int_GW_power_spectra_max);
    arb_init(Int_GW_power_spectra_precision);
    
    arb_init(Int_GW_power_spectra_x_min); // 功率谱积矩形分用
    arb_init(Int_GW_power_spectra_x_max);
    arb_init(Int_GW_power_spectra_x_precision);
    arb_init(Int_GW_power_spectra_y_min);
    arb_init(Int_GW_power_spectra_y_max);
    arb_init(Int_GW_power_spectra_y_precision);
    
}
