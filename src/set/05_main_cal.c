#include "header/05_main_cal.h"

void Set_main_cal(char* comd_argv, slong prec) // comd_argv 为命令行传递参数
{
    //曲率扰动 ζ 相关设定
    // ζ 扰动类型 gaussian_type/exponential_tail_type/up_step_type/power_expansion_type
    Zeta_type=gaussian_type; //需要其它参数之前，如r_m的求解
    
    
    //非高斯相关参数设定，需在求 r_max 和 Mu_2_th 之前设定，其与两都都有关
    //exponential_tail_type
    arb_set_str(Exponential_tail_beta, "-5", prec); //exponential_tail_type 中β的取值，文献中通常取为β=3
    
    /*
    //对应指数尾巴的展开系数 --> (-1)^(n-1)*1/n*β^(n-1)*x^n
    arb_set(Power_expansion_f,Exponential_tail_beta);
    arb_neg(Power_expansion_f,Power_expansion_f);
    arb_div_ui(Power_expansion_f,Power_expansion_f,2,prec);
    arb_pow_ui(Power_expansion_g,Exponential_tail_beta,2,prec);
    arb_div_ui(Power_expansion_g,Power_expansion_g,3,prec);
    arb_pow_ui(Power_expansion_four,Exponential_tail_beta,3,prec);
    arb_neg(Power_expansion_four,Power_expansion_four);
    arb_div_ui(Power_expansion_four,Power_expansion_four,4,prec);
    arb_pow_ui(Power_expansion_five,Exponential_tail_beta,4,prec);
    arb_div_ui(Power_expansion_five,Power_expansion_five,5,prec);
    arb_pow_ui(Power_expansion_six,Exponential_tail_beta,5,prec);
    arb_neg(Power_expansion_six,Power_expansion_six);
    arb_div_ui(Power_expansion_six,Power_expansion_six,6,prec);
    */
    
    //up_step_type
    arb_set(Up_step_h,Upward_step_spectra_h);
    //arb_set_str(Up_step_h, "-2.5", prec); //up_step_type 其中Up_step_h取值为负
    //arb_set_str(Up_step_h, comd_argv, prec); //从命令行读取参数
    
    
    //power_expansion_type 最高可展开到 6 阶
    //通过 h 得到 f=5/12 * |h|
    arb_mul_ui(Power_expansion_f,Up_step_h,5,prec);
    arb_div_ui(Power_expansion_f,Power_expansion_f,12,prec);
    
    //arb_set_str(Power_expansion_f, "-2.5", prec); //power-series expansion 二次项 f_NL -> A
    arb_set_str(Power_expansion_g, "0", prec); //power-series expansion 三次项 g_NL -> B
    arb_set_str(Power_expansion_four, "0", prec); //power-series expansion 四次项 four -> C
    arb_set_str(Power_expansion_five, "0", prec); //power-series expansion 五次项 five -> D
    arb_set_str(Power_expansion_six, "0", prec); //power-series expansion 六次项 six -> E
    //A=B=C=D=E=0 回到高斯的情况
    //printf("command intput arg 1: %s\n",argv[1]);
    
    
    //考虑所有k模式，用δ谱求连续谱，计算∫β(k_◦)dk_◦用
    arb_mul_ui(Log_normal_mul_sigma,Power_sigma, 5, prec);
    arb_sub(PS_abundance_beta_delta_k_all_min,Ln_K_star,Log_normal_mul_sigma,prec); //动态积分范围
    arb_add(PS_abundance_beta_delta_k_all_max,Ln_K_star,Log_normal_mul_sigma,prec);
    arb_set_str(PS_abundance_beta_delta_k_all_precision,"1E-10",prec);
    PS_abundance_beta_delta_k_all_Int_iterate_min=3;
    PS_abundance_beta_delta_k_all_Int_iterate_max=5;
    
    //考虑所有k模式，计算β(M)用
    arb_set_str(PS_abundance_beta_delta_k_M_precision,"1E-20",prec);
    PS_abundance_beta_delta_k_M_Int_iterate_min=5;
    PS_abundance_beta_delta_k_M_Int_iterate_max=10;
    
    
    //设置求最大值区间，如求 r_m Mu_2_th ,需在求ζ(r)参数 µ 的临界值之前设定
    switch(Power_spectrum_type)
    {
        case lognormal_type :
            
            // PS 计算方差 XX YY XY 积分，log-normal谱共用
            //arb_set_str(Power_sigma,"0.03",prec);
            arb_mul_ui(Log_normal_mul_sigma,Power_sigma, 6, prec);
            arb_sub(PS_Int_variance_min,Ln_K_star,Log_normal_mul_sigma,prec); //动态积分范围
            arb_add(PS_Int_variance_max,Ln_K_star,Log_normal_mul_sigma,prec);
            arb_set_str(PS_Int_variance_precision,"1E-20",prec);
            
            switch(Zeta_type)
            {
                case gaussian_type :
                    
                    //arb_set_str(Int_r_min,"1E-14",prec);
                    //arb_set_str(Int_r_max,"9.9E-13",prec);
                    //根据δ谱下r*k为定值，采用动态r_m区间
                    arb_set_str(R_K_to_r_m,"1",prec);
                    arb_div(Int_r_min,R_K_to_r_m,K_star,prec);
                    arb_div_ui(Int_r_min,Int_r_min,10,prec);
                    arb_mul_ui(Int_r_max,Int_r_min,60,prec);
                    Root_r_num=50;
                    arb_set_str(Int_r_precision,"1E-20",prec);
                    
                    arb_set_str(Int_mu_min,"0.3",prec);
                    arb_set_str(Int_mu_max,"1.6",prec);
                    Root_mu_num=25;
                    arb_set_str(Int_mu_precision,"1E-6",prec);
                    
                    C_m_average_iterate_min=3; //求 C_m_average 不好求，迭代次数需单独设置
                    C_m_average_iterate_max=10;
                    arb_set_str(C_m_average_precision,"1E-6",prec);
                    
                    // M -> μ 求根用
                    arb_set_str(Root_M_to_mu_min,"0.1",prec); //Root_M_to_mu_min 最小应该是 Mu_2_th
                    arb_set_str(Root_M_to_mu_max,"2.5",prec);
                    Root_M_to_mu_num=20;
                    arb_set_str(Root_M_to_mu_precision,"1E-6",prec);
                    
                    
                    arb_set_str(Int_n_pk_k_3_min,"0.05",prec); // n_pk(mu_2,k_3) 中 k_3 的积分区间
                    arb_set_str(Int_n_pk_k_3_max,"1.6",prec);
                    arb_set_str(Int_n_pk_k_3_precision,"1E-7",prec);
                    
                    // PS 计算相关 
                    arb_set_str(PS_Int_P_C_l_min,"-1.5",prec); // PS 计算 C_ℓ 的概率密度分布 P(C_l)
                    arb_set_str(PS_Int_P_C_l_max,"1.5",prec);
                    arb_set_str(PS_Int_P_C_l_precision,"1E-40",prec);
                    
                    arb_set_str(PS_abundance_f_all_precision,"1E-10",prec); //PS 最终占比f积分的精度
                    
                    break;
                    
                case exponential_tail_type :
                    
                    //根据δ谱下r*k为定值，采用动态r_m区间
                    arb_set_str(R_K_to_r_m,"2",prec);
                    arb_div(Int_r_min,R_K_to_r_m,K_star,prec);
                    arb_div_ui(Int_r_min,Int_r_min,15,prec);
                    arb_mul_ui(Int_r_max,Int_r_min,80,prec);
                    Root_r_num=80;
                    arb_set_str(Int_r_precision,"1E-20",prec);
                    
                    arb_set_str(Int_mu_min,"0.1",prec);
                    arb_set_str(Int_mu_max,"1.3",prec);
                    Root_mu_num=40;
                    arb_set_str(Int_mu_precision,"1E-6",prec);
                    
                    C_m_average_iterate_min=3; //求 C_m_average 不好求，迭代次数需单独设置
                    C_m_average_iterate_max=7;
                    arb_set_str(C_m_average_precision,"1E-6",prec);
                    
                    // M -> μ 求根用
                    arb_set_str(Root_M_to_mu_min,"0.1",prec); //Root_M_to_mu_min 最小应该是 Mu_2_th
                    arb_set_str(Root_M_to_mu_max,"2.5",prec);
                    Root_M_to_mu_num=20;
                    arb_set_str(Root_M_to_mu_precision,"1E-6",prec);
                    
                    
                    arb_set_str(Int_n_pk_k_3_min,"0.05",prec); // n_pk(mu_2,k_3) 中 k_3 的积分区间
                    arb_set_str(Int_n_pk_k_3_max,"1.6",prec);
                    arb_set_str(Int_n_pk_k_3_precision,"1E-7",prec);
                    
                    // PS 计算相关
                    if( arb_is_positive(Exponential_tail_beta) ) //β>0，定义域 ζ_G<1/β
                    {
                        arb_set_str(PS_Int_P_C_l_min,"-1.5",prec); // PS 计算 C_ℓ 的概率密度分布 P(C_l)
                        arb_set(PS_Int_P_C_l_max,Exponential_tail_beta); //ζ_G<1/β
                        arb_inv(PS_Int_P_C_l_max,PS_Int_P_C_l_max,prec); 
                        arb_set_str(PS_Int_P_C_l_precision,"1E-40",prec);
                    }else //β<0，0，定义域 ζ_G>1/β
                    {
                        arb_set(PS_Int_P_C_l_min,Exponential_tail_beta); //ζ_G>1/β
                        arb_inv(PS_Int_P_C_l_min,PS_Int_P_C_l_min,prec); 
                        arb_set_str(PS_Int_P_C_l_max,"1.5",prec); // PS 计算 C_ℓ 的概率密度分布 P(C_l)
                        arb_set_str(PS_Int_P_C_l_precision,"1E-40",prec);
                    }
                    
                    
                    arb_set_str(PS_abundance_f_all_precision,"1E-10",prec); //PS 最终占比f积分的精度
                    
                    break;
                    
                case power_expansion_type :
                    
                    //根据δ谱下r*k为定值，采用动态r_m区间
                    arb_set_str(R_K_to_r_m,"1",prec);
                    arb_div(Int_r_min,R_K_to_r_m,K_star,prec);
                    arb_div_ui(Int_r_min,Int_r_min,10,prec);
                    arb_mul_ui(Int_r_max,Int_r_min,60,prec);
                    Root_r_num=60;
                    arb_set_str(Int_r_precision,"1E-20",prec);
                    
                    
                    arb_set_str(Int_mu_min,"0.2",prec);
                    arb_set_str(Int_mu_max,"1.3",prec);
                    Root_mu_num=30;
                    arb_set_str(Int_mu_precision,"1E-6",prec);
                    
                    C_m_average_iterate_min=3; //求 C_m_average 不好求，迭代次数需单独设置
                    C_m_average_iterate_max=5;
                    arb_set_str(C_m_average_precision,"1E-6",prec);
                    
                    // M -> μ 求根用
                    arb_set_str(Root_M_to_mu_min,"0.1",prec); //Root_M_to_mu_min 最小应该是 Mu_2_th
                    arb_set_str(Root_M_to_mu_max,"0.7",prec);
                    Root_M_to_mu_num=12;
                    arb_set_str(Root_M_to_mu_precision,"1E-6",prec);
                    
                    
                    arb_set_str(Int_n_pk_k_3_min,"0.05",prec); // n_pk(mu_2,k_3) 中 k_3 的积分区间
                    arb_set_str(Int_n_pk_k_3_max,"1.6",prec);
                    arb_set_str(Int_n_pk_k_3_precision,"1E-7",prec);
                    
                    // PS 计算相关
                    arb_set_str(PS_Int_P_C_l_min,"-1.5",prec); // PS 计算 C_ℓ 的概率密度分布 P(C_l)
                    arb_set_str(PS_Int_P_C_l_max,"1.5",prec);
                    arb_set_str(PS_Int_P_C_l_precision,"1E-50",prec);
                    
                    arb_set_str(PS_abundance_f_all_precision,"1E-10",prec); //PS 最终占比f积分的精度
                    
                    break;
                case up_step_type :
                    
                    //根据δ谱下r*k为定值，采用动态r_m区间
                    arb_set_str(R_K_to_r_m,"9.1",prec);
                    arb_div(Int_r_min,R_K_to_r_m,K_star,prec);
                    arb_div_ui(Int_r_min,Int_r_min,15,prec);
                    arb_mul_ui(Int_r_max,Int_r_min,70,prec);
                    Root_r_num=70;
                    arb_set_str(Int_r_precision,"1E-20",prec);
                    
                    arb_set_str(Int_mu_min,"0.1",prec);
                    arb_set_str(Int_mu_max,"1.5",prec);
                    Root_mu_num=60;
                    arb_set_str(Int_mu_precision,"1E-6",prec);
                    
                    C_m_average_iterate_min=3; //求 C_m_average 不好求，迭代次数需单独设置
                    C_m_average_iterate_max=5;
                    arb_set_str(C_m_average_precision,"1E-6",prec);
                    
                    // M -> μ 求根用
                    arb_set_str(Root_M_to_mu_min,"0.1",prec); //Root_M_to_mu_min 最小应该是 Mu_2_th
                    arb_set_str(Root_M_to_mu_max,"0.7",prec);
                    Root_M_to_mu_num=12;
                    arb_set_str(Root_M_to_mu_precision,"1E-6",prec);
                    
                    
                    arb_set_str(Int_n_pk_k_3_min,"0.05",prec); // n_pk(mu_2,k_3) 中 k_3 的积分区间
                    arb_set_str(Int_n_pk_k_3_max,"1.6",prec);
                    arb_set_str(Int_n_pk_k_3_precision,"1E-7",prec);
                    
                    // PS 计算相关
                    //注意，ζ<2/h，ζ_G<1/h 两者的取值范围不一样
                    arb_set_str(PS_Int_P_C_l_min,"-1.5",prec); // PS 计算 C_ℓ 的概率密度分布 P(C_l)
                    arb_abs(PS_Int_P_C_l_max,Up_step_h);
                    arb_inv(PS_Int_P_C_l_max,PS_Int_P_C_l_max,prec); // ζ_G<1/h
                    arb_set_str(PS_Int_P_C_l_precision,"1E-40",prec);
                    
                    arb_set_str(PS_abundance_f_all_precision,"1E-10",prec); //PS 最终占比f积分的精度
                    
                    break;
                    
                default:
                    printf("main.c Power_spectrum_type->lognormal_type-->zeta_type 有误\n");
                    exit(1);
            }
            
            break;
            
        case delta_type :
            switch(Zeta_type) 
            {
                case gaussian_type :
                    //arb_set_str(Int_r_min,"1E-14",prec);
                    //arb_set_str(Int_r_max,"9.9E-13",prec);
                    //根据δ谱下r*k为定值，采用动态r_m区间
                    arb_set_str(R_K_to_r_m,"1",prec);
                    arb_div(Int_r_min,R_K_to_r_m,K_star,prec);
                    arb_div_ui(Int_r_min,Int_r_min,10,prec);
                    arb_mul_ui(Int_r_max,Int_r_min,60,prec);
                    Root_r_num=60;
                    arb_set_str(Int_r_precision,"1E-23",prec);
                    
                    arb_set_str(Int_mu_min,"0.2",prec);
                    arb_set_str(Int_mu_max,"1.2",prec);
                    Root_mu_num=35;
                    arb_set_str(Int_mu_precision,"1E-9",prec);
                    
                    C_m_average_iterate_min=5; //求 C_m_average 不好求，迭代次数需单独设置
                    C_m_average_iterate_max=10;
                    arb_set_str(C_m_average_precision,"1E-9",prec);
                    
                    
                    // M -> μ 求根用
                    arb_set_str(Root_M_to_mu_min,"0.1",prec); //Root_M_to_mu_min 最小应该是 Mu_2_th
                    arb_set_str(Root_M_to_mu_max,"2.0",prec);
                    Root_M_to_mu_num=20;
                    arb_set_str(Root_M_to_mu_precision,"1E-8",prec);
                    
                    
                    arb_set_str(Int_n_pk_k_3_min,"0.05",prec); // n_pk(mu_2,k_3) 中 k_3 的积分区间
                    arb_set_str(Int_n_pk_k_3_max,"1.6",prec);
                    arb_set_str(Int_n_pk_k_3_precision,"1E-8",prec);
                    
                    // PS 计算相关
                    arb_set_str(PS_abundance_f_all_precision,"1E-10",prec); //PS 最终占比f积分的精度
                    
                    break;
                    
                case power_expansion_type :
                    //arb_set_str(Int_r_min,"5E-15",prec);
                    //arb_set_str(Int_r_max,"9.9E-13",prec);
                    //根据δ谱下r*k为定值，采用动态r_m区间
                    arb_set_str(R_K_to_r_m,"1",prec);
                    arb_div(Int_r_min,R_K_to_r_m,K_star,prec);
                    arb_div_ui(Int_r_min,Int_r_min,10,prec);
                    arb_mul_ui(Int_r_max,Int_r_min,60,prec);
                    Root_r_num=100;
                    arb_set_str(Int_r_precision,"1E-22",prec);
                    
                    
                    arb_set_str(Int_mu_min,"0.2",prec);
                    arb_set_str(Int_mu_max,"1.3",prec);
                    Root_mu_num=40;
                    arb_set_str(Int_mu_precision,"1E-9",prec);
                    
                    
                    C_m_average_iterate_min=5; //求 C_m_average 不好求，迭代次数需单独设置
                    C_m_average_iterate_max=10;
                    arb_set_str(C_m_average_precision,"1E-8",prec);
                    
                    // M -> μ 求根用
                    arb_set_str(Root_M_to_mu_min,"0.1",prec); //Root_M_to_mu_min 最小应该是 Mu_2_th
                    arb_set_str(Root_M_to_mu_max,"2.5",prec);
                    Root_M_to_mu_num=20;
                    arb_set_str(Root_M_to_mu_precision,"1E-7",prec);
                    
                    
                    arb_set_str(Int_n_pk_k_3_min,"0.05",prec); // n_pk(mu_2,k_3) 中 k_3 的积分区间
                    arb_set_str(Int_n_pk_k_3_max,"1.6",prec);
                    arb_set_str(Int_n_pk_k_3_precision,"1E-7",prec);
                    
                    // PS 计算相关
                    arb_set_str(PS_Root_C_l_to_Y_min,"-5",prec); //PS中δ情况，计算P(C_ℓ)需反解Y
                    arb_set_str(PS_Root_C_l_to_Y_max,"5",prec);
                    PS_Root_C_l_to_Y_num=150;
                    arb_set_str(PS_Root_C_l_to_Y_precision,"1E-15",prec);
                    
                    
                    arb_set_str(PS_abundance_f_all_precision,"1E-10",prec); //PS 最终占比f积分的精度
                    
                    break;
                case exponential_tail_type :
                    
                    //arb_set_str(Int_r_min,"1E-14",prec); //小 1.71142578125e-14
                    //arb_set_str(Int_r_max,"9.9E-13",prec);
                    //根据δ谱下r*k为定值，采用动态r_m区间
                    arb_set_str(R_K_to_r_m,"2",prec);
                    arb_div(Int_r_min,R_K_to_r_m,K_star,prec);
                    arb_div_ui(Int_r_min,Int_r_min,15,prec);
                    arb_mul_ui(Int_r_max,Int_r_min,90,prec);
                    Root_r_num=130;
                    arb_set_str(Int_r_precision,"1E-22",prec);
                    
                    arb_set_str(Int_mu_min,"0.1",prec);
                    arb_set_str(Int_mu_max,"1.5",prec);
                    Root_mu_num=50;
                    arb_set_str(Int_mu_precision,"1E-9",prec);
                    
                    C_m_average_iterate_min=3; //求 C_m_average 不好求，迭代次数需单独设置
                    C_m_average_iterate_max=7;
                    arb_set_str(C_m_average_precision,"1E-7",prec);
                    
                    // M -> μ 求根用
                    arb_set_str(Root_M_to_mu_min,"0.1",prec); //Root_M_to_mu_min 最小应该是 Mu_2_th
                    arb_set_str(Root_M_to_mu_max,"2.5",prec);
                    Root_M_to_mu_num=20;
                    arb_set_str(Root_M_to_mu_precision,"1E-7",prec);
                    
                    
                    arb_set_str(Int_n_pk_k_3_min,"0.05",prec); // n_pk(mu_2,k_3) 中 k_3 的积分区间
                    arb_set_str(Int_n_pk_k_3_max,"1.6",prec);
                    arb_set_str(Int_n_pk_k_3_precision,"1E-7",prec);
                    
                    // PS 计算相关
                    arb_set_str(PS_Root_C_l_to_Y_min,"-5",prec); //PS中δ情况，计算P(C_ℓ)需反解Y
                    arb_set_str(PS_Root_C_l_to_Y_max,"5",prec);
                    PS_Root_C_l_to_Y_num=150;
                    arb_set_str(PS_Root_C_l_to_Y_precision,"1E-15",prec);
                    
                    arb_set_str(PS_abundance_f_all_precision,"1E-10",prec); //PS 最终占比f积分的精度
                    
                    break;
                    
                case up_step_type :
                    
                    //arb_set_str(Int_r_min,"1E-14",prec); //小 1.71142578125e-14
                    //arb_set_str(Int_r_max,"9.9E-13",prec);
                    //根据δ谱下r*k为定值，采用动态r_m区间
                    arb_set_str(R_K_to_r_m,"9",prec);
                    arb_div(Int_r_min,R_K_to_r_m,K_star,prec);
                    arb_div_ui(Int_r_min,Int_r_min,15,prec);
                    arb_mul_ui(Int_r_max,Int_r_min,80,prec);
                    Root_r_num=120;
                    arb_set_str(Int_r_precision,"1E-22",prec);
                    
                    arb_set_str(Int_mu_min,"0.1",prec);
                    arb_set_str(Int_mu_max,"1.5",prec);
                    Root_mu_num=45;
                    arb_set_str(Int_mu_precision,"1E-9",prec);
                    
                    C_m_average_iterate_min=3; //求 C_m_average 不好求，迭代次数需单独设置
                    C_m_average_iterate_max=5;
                    arb_set_str(C_m_average_precision,"1E-7",prec);
                    
                    // M -> μ 求根用
                    arb_set_str(Root_M_to_mu_min,"0.1",prec); //Root_M_to_mu_min 最小应该是 Mu_2_th
                    arb_set_str(Root_M_to_mu_max,"0.7",prec);
                    Root_M_to_mu_num=10;
                    arb_set_str(Root_M_to_mu_precision,"1E-7",prec);
                    
                    
                    arb_set_str(Int_n_pk_k_3_min,"0.05",prec); // n_pk(mu_2,k_3) 中 k_3 的积分区间
                    arb_set_str(Int_n_pk_k_3_max,"1.6",prec);
                    arb_set_str(Int_n_pk_k_3_precision,"1E-7",prec);
                    
                    // PS 计算相关
                    arb_set_str(PS_Root_C_l_to_Y_min,"-5",prec); //PS中δ情况，计算P(C_ℓ)需反解Y
                    arb_set_str(PS_Root_C_l_to_Y_max,"5",prec);
                    PS_Root_C_l_to_Y_num=150;
                    arb_set_str(PS_Root_C_l_to_Y_precision,"1E-15",prec);
                    
                    arb_set_str(PS_abundance_f_all_precision,"1E-10",prec); //PS 最终占比f积分的精度
                    
                    break;
                    
                default :
                    printf("main.c Power_spectrum_type->delta_type-->zeta_type 有误\n");
                    exit(1);
            }
                    
            break;
        case power_law_type :
            switch(Zeta_type) 
            {
                case gaussian_type:
                    arb_set_str(Int_r_min,"1E-1",prec);
                    arb_set_str(Int_r_max,"10",prec);
                    Root_r_num=25;
                    arb_set_str(Int_r_precision,"1E-20",prec);
                    
                    arb_set_str(Int_mu_min,"0.4",prec);
                    arb_set_str(Int_mu_max,"1.5",prec);
                    Root_mu_num=15;
                    arb_set_str(Int_mu_precision,"1E-9",prec);
                    
                    C_m_average_iterate_min=3; //求 C_m_average 不好求，迭代次数需单独设置
                    C_m_average_iterate_max=10;
                    arb_set_str(C_m_average_precision,"1E-7",prec);
                    
                    // M -> μ 求根用
                    arb_set_str(Root_M_to_mu_min,"0.1",prec); //Root_M_to_mu_min 最小应该是 Mu_2_th
                    arb_set_str(Root_M_to_mu_max,"2.5",prec);
                    Root_M_to_mu_num=10;
                    arb_set_str(Root_M_to_mu_precision,"1E-7",prec);
                    
                    
                    arb_set_str(Int_n_pk_k_3_min,"0.05",prec); // n_pk(mu_2,k_3) 中 k_3 的积分区间
                    arb_set_str(Int_n_pk_k_3_max,"1.6",prec);
                    arb_set_str(Int_n_pk_k_3_precision,"1E-7",prec);
                    
                    // PS 计算相关 
                    //arb_set_str(Power_sigma,"0.03",prec);
                    arb_sub_ui(PS_Int_variance_min,Ln_K_star,15,prec); // PS 计算方差 XX YY XY 积分
                    arb_add_ui(PS_Int_variance_max,Ln_K_star,15,prec);
                    arb_set_str(PS_Int_variance_precision,"1E-20",prec);
                    
                    arb_set_str(PS_Int_P_C_l_min,"-1.5",prec); // PS 计算 C_ℓ 的概率密度分布 P(C_l)
                    arb_set_str(PS_Int_P_C_l_max,"1.5",prec);
                    arb_set_str(PS_Int_P_C_l_precision,"1E-40",prec);
                    
                    arb_set_str(PS_abundance_f_all_precision,"1E-10",prec); //PS 最终占比f积分的精度
                    
                    break;
                default :
                    printf("main.c Power_spectrum_type->power_law_type-->zeta_type 有误\n");
                    exit(1);
            }
            break;
        case broken_power_law_type :
            // PS 计算方差 XX YY XY 积分， broken_power_law_type 谱共用
            //在左侧大约以 α 斜率上升，右侧大约经 β 斜率下降，跟据 α/β 值动态取值
            //e^(-n)*k^α=K_star^α => lnk = ln(k_star) -n/α
            arb_ui_div(PS_Int_variance_min,28,BPL_alpha,prec); //左边，n=28
            arb_sub(PS_Int_variance_min,Ln_K_star,PS_Int_variance_min,prec);
            
            arb_ui_div(PS_Int_variance_max,28,BPL_beta,prec); //右边，n=28
            arb_add(PS_Int_variance_max,Ln_K_star,PS_Int_variance_max,prec);
            
            arb_set_str(PS_Int_variance_precision,"1E-19",prec);
            
            switch(Zeta_type)
            {
                case gaussian_type :
                    
                    //根据δ谱下r*k为定值，采用动态r_m区间
                    arb_set_str(R_K_to_r_m,"3",prec);
                    arb_div(Int_r_min,R_K_to_r_m,K_star,prec);
                    arb_div_ui(Int_r_min,Int_r_min,15,prec);
                    arb_mul_ui(Int_r_max,Int_r_min,80,prec);
                    Root_r_num=120;
                    arb_set_str(Int_r_precision,"1E-20",prec);
                    
                    arb_set_str(Int_mu_min,"0.3",prec);
                    arb_set_str(Int_mu_max,"1.2",prec);
                    Root_mu_num=15;
                    arb_set_str(Int_mu_precision,"1E-6",prec);
                    
                    C_m_average_iterate_min=3; //求 C_m_average 不好求，迭代次数需单独设置
                    C_m_average_iterate_max=10;
                    arb_set_str(C_m_average_precision,"1E-6",prec);
                    
                    // M -> μ 求根用
                    arb_set_str(Root_M_to_mu_min,"0.1",prec); //Root_M_to_mu_min 最小应该是 Mu_2_th
                    arb_set_str(Root_M_to_mu_max,"2.5",prec);
                    Root_M_to_mu_num=20;
                    arb_set_str(Root_M_to_mu_precision,"1E-6",prec);
                    
                    
                    arb_set_str(Int_n_pk_k_3_min,"0.05",prec); // n_pk(mu_2,k_3) 中 k_3 的积分区间
                    arb_set_str(Int_n_pk_k_3_max,"1.6",prec);
                    arb_set_str(Int_n_pk_k_3_precision,"1E-7",prec);
                    
                    // PS 计算相关
                    arb_set_str(PS_Int_P_C_l_min,"-1.5",prec); // PS 计算 C_ℓ 的概率密度分布 P(C_l)
                    arb_set_str(PS_Int_P_C_l_max,"1.5",prec);
                    arb_set_str(PS_Int_P_C_l_precision,"1E-40",prec);
                    
                    arb_set_str(PS_abundance_f_all_precision,"1E-10",prec); //PS 最终占比f积分的精度
                    
                    break;
                case up_step_type :
                    
                    //根据δ谱下r*k为定值，采用动态r_m区间
                    arb_set_str(R_K_to_r_m,"9.1",prec);
                    arb_div(Int_r_min,R_K_to_r_m,K_star,prec);
                    arb_div_ui(Int_r_min,Int_r_min,15,prec);
                    arb_mul_ui(Int_r_max,Int_r_min,70,prec);
                    Root_r_num=70;
                    arb_set_str(Int_r_precision,"1E-20",prec);
                    
                    arb_set_str(Int_mu_min,"0.1",prec);
                    arb_set_str(Int_mu_max,"1.5",prec);
                    Root_mu_num=60;
                    arb_set_str(Int_mu_precision,"1E-6",prec);
                    
                    C_m_average_iterate_min=3; //求 C_m_average 不好求，迭代次数需单独设置
                    C_m_average_iterate_max=5;
                    arb_set_str(C_m_average_precision,"1E-6",prec);
                    
                    // M -> μ 求根用
                    arb_set_str(Root_M_to_mu_min,"0.1",prec); //Root_M_to_mu_min 最小应该是 Mu_2_th
                    arb_set_str(Root_M_to_mu_max,"0.7",prec);
                    Root_M_to_mu_num=12;
                    arb_set_str(Root_M_to_mu_precision,"1E-6",prec);
                    
                    
                    arb_set_str(Int_n_pk_k_3_min,"0.05",prec); // n_pk(mu_2,k_3) 中 k_3 的积分区间
                    arb_set_str(Int_n_pk_k_3_max,"1.6",prec);
                    arb_set_str(Int_n_pk_k_3_precision,"1E-7",prec);
                    
                    // PS 计算相关
                    //注意，ζ<2/h，ζ_G<1/h 两者的取值范围不一样
                    arb_set_str(PS_Int_P_C_l_min,"-1.5",prec); // PS 计算 C_ℓ 的概率密度分布 P(C_l)
                    arb_abs(PS_Int_P_C_l_max,Up_step_h);
                    arb_inv(PS_Int_P_C_l_max,PS_Int_P_C_l_max,prec); // ζ_G<1/h
                    arb_set_str(PS_Int_P_C_l_precision,"1E-35",prec);
                    
                    arb_set_str(PS_abundance_f_all_precision,"1E-10",prec); //PS 最终占比f积分的精度
                    
                    break;
                default :
                    printf("main.c Power_spectrum_type->broken_power_law_type-->zeta_type 有误\n");
                    exit(1);
            }
            
            break;
        case box_type :
            switch(Zeta_type)
            {
                case gaussian_type:
                    arb_set_str(Int_r_min,"1E-1",prec);
                    arb_set_str(Int_r_max,"9.9",prec);
                    Root_r_num=50;
                    arb_set_str(Int_r_precision,"1E-20",prec);
                    
                    arb_set_str(Int_mu_min,"0.55",prec);
                    arb_set_str(Int_mu_max,"1.5",prec);
                    Root_mu_num=10;
                    arb_set_str(Int_mu_precision,"1E-9",prec);
                    
                    C_m_average_iterate_min=3; //求 C_m_average 不好求，迭代次数需单独设置
                    C_m_average_iterate_max=10;
                    arb_set_str(C_m_average_precision,"1E-7",prec);
                    
                    // M -> μ 求根用
                    arb_set_str(Root_M_to_mu_min,"0.1",prec); //Root_M_to_mu_min 最小应该是 Mu_2_th
                    arb_set_str(Root_M_to_mu_max,"2.5",prec);
                    Root_M_to_mu_num=10;
                    arb_set_str(Root_M_to_mu_precision,"1E-7",prec);
                    
                    
                    arb_set_str(Int_n_pk_k_3_min,"0.05",prec); // n_pk(mu_2,k_3) 中 k_3 的积分区间
                    arb_set_str(Int_n_pk_k_3_max,"1.6",prec);
                    arb_set_str(Int_n_pk_k_3_precision,"1E-7",prec);
                    
                    // PS 计算相关 
                    arb_set(PS_Int_variance_min,Int_sigma_n_min); // PS 计算方差 XX YY XY 积分
                    arb_set(PS_Int_variance_max,Int_sigma_n_max);
                    arb_set_str(PS_Int_variance_precision,"1E-20",prec);
                    
                    arb_set_str(PS_Int_P_C_l_min,"-1.5",prec); // PS 计算 C_ℓ 的概率密度分布 P(C_l)
                    arb_set_str(PS_Int_P_C_l_max,"1.5",prec);
                    arb_set_str(PS_Int_P_C_l_precision,"1E-40",prec);
                    
                    arb_set_str(PS_abundance_f_all_precision,"1E-10",prec); //PS 最终占比f积分的精度
                    
                    break;
                default :
                    printf("main.c Power_spectrum_type->box_type-->zeta_type 有误\n");
                    exit(1);
            }
            break;
        case link_cmb_type :
            switch(Zeta_type) 
            {
                case gaussian_type:
                    arb_set_str(Int_r_min,"1E-16",prec);
                    arb_set_str(Int_r_max,"9.9E-15",prec);
                    Root_r_num=100;
                    arb_set_str(Int_r_precision,"1E-20",prec);
                    
                    arb_set_str(Int_mu_min,"0.1",prec);
                    arb_set_str(Int_mu_max,"1.2",prec);
                    Root_mu_num=10;
                    arb_set_str(Int_mu_precision,"1E-9",prec);
                    
                    C_m_average_iterate_min=3; //求 C_m_average 不好求，迭代次数需单独设置
                    C_m_average_iterate_max=10;
                    arb_set_str(C_m_average_precision,"1E-7",prec);
                    
                    // M -> μ 求根用
                    arb_set_str(Root_M_to_mu_min,"0.1",prec); //Root_M_to_mu_min 最小应该是 Mu_2_th
                    arb_set_str(Root_M_to_mu_max,"2.5",prec);
                    Root_M_to_mu_num=10;
                    arb_set_str(Root_M_to_mu_precision,"1E-7",prec);
                    
                    
                    arb_set_str(Int_n_pk_k_3_min,"0.05",prec); // n_pk(mu_2,k_3) 中 k_3 的积分区间
                    arb_set_str(Int_n_pk_k_3_max,"1.6",prec);
                    arb_set_str(Int_n_pk_k_3_precision,"1E-7",prec);
                    
                    // PS 计算相关 
                    arb_set(PS_Int_variance_min,Int_sigma_n_min); // PS 计算方差 XX YY XY 积分
                    arb_set(PS_Int_variance_max,Int_sigma_n_max);
                    arb_set_str(PS_Int_variance_precision,"1E-20",prec);
                    
                    arb_set_str(PS_Int_P_C_l_min,"-1.5",prec); // PS 计算 C_ℓ 的概率密度分布 P(C_l)
                    arb_set_str(PS_Int_P_C_l_max,"1.5",prec);
                    arb_set_str(PS_Int_P_C_l_precision,"1E-40",prec);
                    
                    arb_set_str(PS_abundance_f_all_precision,"1E-10",prec); //PS 最终占比f积分的精度
                    
                    break;
                default :
                    printf("main.c Power_spectrum_type->box_type-->zeta_type 有误\n");
                    exit(1);
            }
            break;
        case upward_step_spectra_type :
            //并不做PBHs相关的具体计算，仅计算SIGWs，随意设
            // PS 计算方差 XX YY XY 积分， upward_step_spectra_type 谱共用
            arb_set_str(PS_Int_variance_min, "0", prec);
            arb_set_str(PS_Int_variance_max, "0", prec);
            arb_set_str(PS_Int_variance_precision,"0",prec);
            
            break;
        default :
            printf(" main.c Power_spectrum_type 有误\n");
            exit(1);
    }
    
    
    //求解 Mu_2_th 相关设定
    
    //变量 Mu_2，ζ(r) 的参数 µ
    arb_one(Mu_2); //ζ(r) 的参数 µ，默认为Mu_2=1
    
    //变量 K_3_square，ζ(r) 的参数 k
    arb_one(K_3_square); //默认为K_3_square=1
    
    //此设定应在计算 Mu_2_th 之前，有此设定有影响 Mu_2_th 的值
    FIT_FUNC_IF=false; //在后需的计算中，若拟合完成，是否开启拟合 true/false
    Relative_Mass=true; //计算黑洞的质量分布时，是否使用相对质量来进行表示和计算
                        //即： β(M)-->β(M/M_H) ， f(M)-->f(M/M_H)
                        // 非相对质量的计算还有点小问题
    SIMPLIFY=true; //是否启用简化版本的计算
    //简化启用时，K_3_square设为其平均值，此时ζ^G(r)的表达式大为化简，只剩下一个随机变量 μ_2
    
    if( SIMPLIFY )
    {
        arb_set(K_3_square,Gamma_3); // 设 K_3_square=Gamma_3 时，
        //极大化简 ζ(r) 的表达式 ζ(r)=Mu_2 * ψ_1(r)
        //并且，此时 μ_2th 也不会变动 (K_3_square 改变会导致 μ_2th 改变)
    }
    
    
    //在视界进入时，视界质量 M_H，形成黑洞质量 M，两者间的关系可近似看作 scaling law 形式
    //临界坍缩： M=K*(C-C_th)^γ * M_H
    //arb_one(Mass_K); //这里设 K=1
    arb_set_str(Mass_K,"4",prec); //这里设 K=4
    
    //arb_set_str(Mass_gamma,"0.36",prec); // γ ≃ 0.36
    arb_set_str(Mass_gamma,"0.357",prec); // γ ≃ 0.357
    
    //β 到 f 的转换系数，可以采用不同方法 beta_f_general_I, beta_f_general_II, beta_f_myself
    beta_to_f_type=beta_f_general_I;
    
    
    //
    //诱导引力波相关设定
    //
    //诱导引力波相关计算方法设定
    // 单就计算速度而言， 功率谱（两种结果不一样） Espinosa_01，能量密度谱(两种结果相同) Kohri_02
    GW_induced_method=Li_gauss; // Li_gauss, Kohri_02, Espinosa_01
    
    
    //诱导引力波积分辐助函数的积分设定
    arb_zero(Int_GW_I_func_min); //积分下限，为 1/k or 0 or 1，设为0可以在很大程度上加快计算
    arb_pos_inf(Int_GW_I_func_max); //积分上限 +∞
    //arb_set_str(Int_GW_I_func_max,"1E4",prec);
    arb_set_str(Int_GW_I_func_precision,"1E-15",prec);
    Int_GW_I_func_iterate_min=7;
    Int_GW_I_func_iterate_max=15;
    
    //诱导引力波功率谱的积分设定
    arb_set_str(Int_GW_power_spectra_min,"1E-15",prec); //为f(x,y)二元积分，这里设定的是x积分的范围，x的下界为0
    arb_set_str(Int_GW_power_spectra_max,"6E3",prec); //上界理论上为 +∞ 
    arb_set_str(Int_GW_power_spectra_precision,"1E-15",prec);
    Int_GW_power_spectra_iterate_min=8;
    Int_GW_power_spectra_iterate_max=13;
    
    //前面诱导引力波的参数设定为是Kohri_02, Espinosa_01两种方法准备的
    //现建议采用 Li_gauss 的方法，其为矩形积分，速度有很大提高，且方便算非高斯性
    //虽然不采用 adaptive 且用 double_exponential 的积分方式，对于一些点而言，速度非常快
    //但考虑到这种情况，对于某些点（如log-normal峰值处的点）非常难算，综合考虑
    //采用 adaptive 方式且此时的积分采用 gauss_kronrod_iterate 
    //且在计算时，不管前面的积分方式是什么，需要临时采用 gauss_kronrod_iterate
    
    // 功率谱积矩形分用
    arb_set_str(Int_GW_power_spectra_x_min,"0",prec); //对x积分 [0, +∞]
    arb_set_str(Int_GW_power_spectra_x_max,"3E3",prec);
    arb_set_str(Int_GW_power_spectra_x_precision,"1E-15",prec);
    
    arb_set_str(Int_GW_power_spectra_y_min,"-1",prec); //对y积分 [-1, 1]
    arb_set_str(Int_GW_power_spectra_y_max,"1",prec);
    arb_set_str(Int_GW_power_spectra_y_precision,"1E-15",prec);
    
    Int_GW_power_spectra_rectangle_adaptive=true;
    
    if( Int_GW_power_spectra_rectangle_adaptive )
    {
        Int_GW_power_spectra_iterate_x_min=6; //x精度
        Int_GW_power_spectra_iterate_x_max=90;
        
        Int_GW_power_spectra_iterate_y_min=4; //y精度
        Int_GW_power_spectra_iterate_y_max=90;
    }else
    {
        Int_GW_power_spectra_iterate_x_min=6; //x精度
        Int_GW_power_spectra_iterate_x_max=14;
        
        Int_GW_power_spectra_iterate_y_min=5; //y精度
        Int_GW_power_spectra_iterate_y_max=12;
    }
    
    
    //诱导引力波中的非高斯项计算
    //没做线程限制，会自动进行多线程运算
    GW_dim_integral_res_print=false; //积分结果打印
    GW_dim_8_MINEVAL=250000; //最小计算次数
    GW_dim_8_MAXEVAL=30000000; //最大计算次数
    GW_dim_8_NSTART=10000; //每次迭代计算次数
    GW_dim_8_NINCREASE=5000; //迭代计算增加次数
    GW_dim_8_EPSREL=1e-3; //相对精度
    GW_dim_8_EPSABS=1e-15; //绝对精度
    GW_dim_8_t_upper="1.5E2"; //积分区间 t 变量上限，理论上为 +∞
    
    GW_dim_6_MINEVAL=99000;
    GW_dim_6_MAXEVAL=10000000;
    GW_dim_6_NSTART=8000;
    GW_dim_6_NINCREASE=4000;
    GW_dim_6_EPSREL=1e-4;
    GW_dim_6_EPSABS=1e-15;
    GW_dim_6_t_upper="3E2";
    
    GW_dim_5_MINEVAL=60000;
    GW_dim_5_MAXEVAL=6000000;
    GW_dim_5_NSTART=5000;
    GW_dim_5_NINCREASE=2500;
    GW_dim_5_EPSREL=1e-5;
    GW_dim_5_EPSABS=1e-15;
    GW_dim_5_t_upper="4E2";
    
    GW_dim_4_MINEVAL=30000;
    GW_dim_4_MAXEVAL=4000000;
    GW_dim_4_NSTART=2000;
    GW_dim_4_NINCREASE=1500;
    GW_dim_4_EPSREL=1e-6;
    GW_dim_4_EPSABS=1e-15;
    GW_dim_4_t_upper="1E3";
    
    GW_dim_2_MINEVAL=6000;
    GW_dim_2_MAXEVAL=600000;
    GW_dim_2_NSTART=1000;
    GW_dim_2_NINCREASE=500;
    GW_dim_2_EPSREL=1e-8;
    GW_dim_2_EPSABS=1e-15;
    GW_dim_2_t_upper="6E3";
}
