#ifndef __PBHS_PHY_CONSTANT__   /* Include guard */
#define __PBHS_PHY_CONSTANT__

#include "new_type.h" //自定义的新类型
#include <arb.h> //高精度实数运算
#include <stdbool.h> // bool 类型
#include <time.h> //用于计算程序耗时
#include <limits.h> //用于获取系统处理相关信息

//使用 extern 声明为外部变量，达到全局变量目的
//还需要在.c模快中赋值

//功率谱相关
extern POWER Power_spectrum_type; //功率谱类型
extern POWER Continuum_spectrum_type;
extern arb_t Power_A; //功率谱振幅
extern arb_t Power_sigma; //功率谱掌宽，例如 log-normal 的
extern arb_t Log_normal_mul_sigma; //功率谱掌宽的倍数
extern arb_t K_star; //功率谱参考尺度
extern arb_t Ln_K_star; //功率谱参考尺度，取对数
extern arb_t Box_K_middle; //Box功率谱 中心点
extern arb_t Box_width; //Box功率谱 宽度
extern arb_t BPL_alpha; //broken power law功率谱使用 α
extern arb_t BPL_beta; //broken power law功率谱使用 β
extern arb_t BPL_gamma; //broken power law功率谱使用 γ
extern arb_t Link_CMB_K_t; //link cmb 功率谱用
extern arb_t Link_CMB_K_m;
extern arb_t Link_CMB_A_s;
extern arb_t Link_CMB_A_b;
extern arb_t Link_CMB_K_star;
extern arb_t Link_CMB_P_m;
extern arb_t Link_CMB_n_s;
extern arb_t Link_CMB_n_b;
extern arb_t Upward_step_spectra_epsilon_S; //upward step模型的原生功率谱, SR-USR-SR
extern arb_t Upward_step_spectra_epsilon_V;
extern arb_t Upward_step_spectra_eta_V;
extern arb_t Upward_step_spectra_g;
extern arb_t Upward_step_spectra_h;
extern arb_t Upward_step_spectra_V_0;
extern arb_t Upward_step_spectra_Delta_V;
extern arb_t Upward_step_spectra_Delta_phi_USR;
extern arb_t Upward_step_spectra_phi_s;
extern arb_t Upward_step_spectra_phi_c;
extern arb_t Upward_step_spectra_k_s;
extern arb_t Upward_step_spectra_k_c;
extern arb_t Upward_step_spectra_ln_k_s;
extern arb_t Upward_step_spectra_ln_k_c;
extern arb_t Upward_step_spectra_tau_c;

//曲率扰动 ζ 相关
extern ZETA Zeta_type; // ζ 类型
extern arb_t PT_mu; // ζ(r) 的参数µ
extern arb_t PT_mu_th; // ζ(r) 的参数 µ 的临界值
extern arb_t PT_mu_max; // ζ(r) 的参数 µ 的最大值
extern arb_t PT_k; // ζ(r) 的参数k
extern arb_t PT_beta_cal_need_exp_2_zeta_m; // 某个 M 对应的 e^(2*ζ_m)
extern arb_t R_MAX; // ζ(r) 取最大值时的 r 值
extern arb_t R_K_to_r_m; //求r_m动态区间用
extern arb_t Q_parameter_th; //ζ(r)取临界值时的q参数
extern arb_t R_m_times_K; //对于δ谱，r_m*k_*为一定值


//计算ζ(r)的辅助变量
extern arb_t Sigma_0_square; //(σ_n)^2
extern arb_t Sigma_1_square;
extern arb_t Sigma_2_square;
extern arb_t Sigma_3_square;
extern arb_t Sigma_4_square;
extern arb_t Gamma_1; //(γ_1)
extern arb_t Gamma_3; //(γ_3)
extern arb_t R_1; //(R_1)
extern arb_t R_3; //(R_3)

extern bool Peak_theory_source_zeta_gradient; //peak theory计算中，统计量是否取ζ的梯度
 
// C(r) 相关
extern ulong C_m_average_iterate_min; //求 C_m_average
extern ulong C_m_average_iterate_max;
extern arb_t C_m_average_precision; 


extern enum PT_MU_TH_enum PT_Mu_th_METHOD; ///求临界值的方法


//在视界进入时，视界质量M_H，形成的黑洞的质量为M，两者间的关系
//与曲率扰动 K 有关，与 γ 也有关
extern arb_t Mass_K; //曲率K
extern arb_t Mass_gamma; //γ
extern bool Critical_Collapse_Effect; //是否考虑临界坍缩效应


//宇宙学基本常数
extern arb_t Omega_DM; //暗物质所占比例
extern arb_t Omega_M; //baryons + dark matter 所占比例
extern arb_t Omega_B; //普通重子物质 baryons
extern arb_t Omega_Lambda; //暗能量 Ω_Λ
extern arb_t Omega_radiation; //辐射所占的比例 radiation density
extern arb_t Omega_gamma; //光子所占的比例 photon density
extern arb_t Hubble_constant_h; //哈勃常数 h
extern arb_t Hubble_constant_H_0; //哈勃常数 H_0
extern arb_t Scale_factor_a; // 尺度因子 a(t)，某些公式的需要
extern arb_t Scale_factor_a_0; //尺度因子 a(t_0)
extern arb_t Equation_Of_State_w; //状态方程 w
extern arb_t K_scale_eq; //在物质辐射相等时刻，对应的共动波数
extern arb_t Z_scale_eq; //在物质辐射相等时刻，对应的红移
extern arb_t A_scale_eq; //在物质辐射相等时刻，对应的尺度因子
extern arb_t T_scale_eq; //在物质辐射相等时刻，对应的温度

extern arb_t effective_g_star_eq; //在物质辐射相等时刻，对应的相对论有效自由度数目
extern arb_t effective_g_star_eq_entropy; //在物质辐射相等时刻，对应的熵有效自由度数目
extern arb_t effective_g_star; //重新进入视界后，对应的相对论有效自由度数目
extern arb_t effective_g_star_entropy; //重新进入视界后，对应的熵有效自由度数目
extern arb_t effective_g_star_current; //当今对应的相对论有效自由度数目
extern arb_t effective_g_star_current_entropy; //当今对应的熵有效自由度数目

//设置求值区间，如求 r_m PT_mu_th
extern arb_t Int_r_min;
extern arb_t Int_r_max;
extern slong Root_r_num;
extern arb_t Int_r_precision;
extern arb_t  Int_mu_min;
extern arb_t  Int_mu_max;
extern slong  Root_mu_num;
extern arb_t Int_mu_precision;

extern arb_t Int_sigma_n_min; //sigma_n
extern arb_t Int_sigma_n_max;
extern arb_t Int_sigma_n_precision;


extern arb_t Int_n_pk_k_min; // n_pk(mu,k) 中 k 的积分区间
extern arb_t Int_n_pk_k_max;
extern arb_t Int_n_pk_k_precision;



extern arb_t Root_M_to_mu_min; // M -> μ 求根用
extern arb_t Root_M_to_mu_max;
extern slong Root_M_to_mu_num;
extern arb_t Root_M_to_mu_precision;

extern ulong Integration_iterate_min; //积分的最小迭代次数
extern ulong Integration_iterate_max; //积分的最大迭代次数


// β 与 f 间的转换系数的选择
extern enum BETA_TO_F_ENUM_TYPE beta_to_f_type;


//非高斯相关参数
extern arb_t Power_expansion_f; //power-series expansion 二次项 f_NL -> A
extern arb_t Power_expansion_g; //power-series expansion 三次项 g_NL -> B
extern arb_t Power_expansion_four; //power-series expansion 四次项 four -> C
extern arb_t Power_expansion_five; //power-series expansion 五次项 five -> D
extern arb_t Power_expansion_six; //power-series expansion 五次项 six -> E
extern arb_t Up_step_h; //up-step 模型中用
extern arb_t Exponential_tail_beta; //exponential_tail 模型中用


//诱导引力波
extern arb_t Int_GW_I_func_min; // I/I_c/I_s 积分用
extern arb_t Int_GW_I_func_max;
extern arb_t Int_GW_I_func_precision;
extern ulong Int_GW_I_func_iterate_min;
extern ulong Int_GW_I_func_iterate_max;

extern arb_t Int_GW_power_spectra_min; // 功率谱积分用
extern arb_t Int_GW_power_spectra_max;
extern arb_t Int_GW_power_spectra_precision;
extern ulong Int_GW_power_spectra_iterate_min;
extern ulong Int_GW_power_spectra_iterate_max;

extern arb_t Int_GW_power_spectra_x_min; // 功率谱积矩形分用
extern arb_t Int_GW_power_spectra_x_max;
extern arb_t Int_GW_power_spectra_x_precision;
extern ulong Int_GW_power_spectra_iterate_x_min;
extern ulong Int_GW_power_spectra_iterate_x_max;
extern arb_t Int_GW_power_spectra_y_min;
extern arb_t Int_GW_power_spectra_y_max;
extern arb_t Int_GW_power_spectra_y_precision;
extern ulong Int_GW_power_spectra_iterate_y_min;
extern ulong Int_GW_power_spectra_iterate_y_max;

extern bool Int_GW_power_spectra_rectangle_adaptive; //矩形二维积分是否采用自适应
extern enum GW_INDUCED_ENUM_TYPE GW_induced_method; //诱导引力波计算方法

extern bool GW_dim_integral_res_print;
extern long int GW_dim_8_MINEVAL; //诱导引力波非高斯性积分参数
extern long int GW_dim_8_MAXEVAL;
extern long int GW_dim_8_NSTART;
extern long int GW_dim_8_NINCREASE;
extern long double GW_dim_8_EPSREL;
extern long double GW_dim_8_EPSABS;
extern char* GW_dim_8_t_upper;

extern long int GW_dim_6_MINEVAL;
extern long int GW_dim_6_MAXEVAL;
extern long int GW_dim_6_NSTART;
extern long int GW_dim_6_NINCREASE;
extern long double GW_dim_6_EPSREL;
extern long double GW_dim_6_EPSABS;
extern char* GW_dim_6_t_upper;

extern long int GW_dim_5_MINEVAL;
extern long int GW_dim_5_MAXEVAL;
extern long int GW_dim_5_NSTART;
extern long int GW_dim_5_NINCREASE;
extern long double GW_dim_5_EPSREL;
extern long double GW_dim_5_EPSABS;
extern char* GW_dim_5_t_upper;

extern long int GW_dim_4_MINEVAL;
extern long int GW_dim_4_MAXEVAL;
extern long int GW_dim_4_NSTART;
extern long int GW_dim_4_NINCREASE;
extern long double GW_dim_4_EPSREL;
extern long double GW_dim_4_EPSABS;
extern char* GW_dim_4_t_upper;

extern long int GW_dim_2_MINEVAL;
extern long int GW_dim_2_MAXEVAL;
extern long int GW_dim_2_NSTART;
extern long int GW_dim_2_NINCREASE;
extern long double GW_dim_2_EPSREL;
extern long double GW_dim_2_EPSABS;
extern char* GW_dim_2_t_upper;

//数学计算
extern arb_t Pi; //常数π
extern arb_t Pi_2; //常数2π (πx2)
extern slong Multithreaded_divide_integration_interval_number; //多线程计算区间分隔数目
extern slong Multithreaded_number; //多线程计算线程数

//extern slong prec; //控制计算精度
//extern slong N_n; //整数n，用于计算辅助

///extern arb_t _TEMP_arb_;//全局临时传递使用 arb_t 
//extern slong _TEMP_slong_;//全局临时传递使用 slong


//函数数据的输入输出，及拟合

extern arb_t Func_output_x_min; //输出范围
extern arb_t Func_output_x_max; 
extern ulong Func_output_number; //输出点的数目
extern char* Func_output_file ; //输出文件名
extern char* Func_output_delimiter; //输出数据间的分格符
extern unsigned int Func_output_skip_header; //输出中添加头开始几行的信息说明
extern unsigned int Func_output_column_x; //可能有多列数据，表明位置
extern unsigned int Func_output_column_y;
extern char * Func_output_fitted_file; //利用输出数据，拟合结果所存储名
extern size_t Fit_entries_number; //拟合读取数据的点数

//函数拟合相关
extern bool FIT_FUNC_IF; //是否开启拟合
extern struct FUNC_FITTED_DATE* Fit_test; //测试用

extern struct FUNC_FITTED_DATE* Fit_psi_1_0; // ψ_1(r) 及其各阶导数 拟合结果
extern struct FUNC_FITTED_DATE* Fit_psi_1_1;
extern struct FUNC_FITTED_DATE* Fit_psi_1_2;
extern struct FUNC_FITTED_DATE* Fit_psi_1_3;
extern struct FUNC_FITTED_DATE* Fit_psi_1_4;

extern struct FUNC_FITTED_DATE* Fit_Laplacian_psi_1_0; // ∆ψ_1(r) 及其各阶导数 拟合结果
extern struct FUNC_FITTED_DATE* Fit_Laplacian_psi_1_1;
extern struct FUNC_FITTED_DATE* Fit_Laplacian_psi_1_2;
extern struct FUNC_FITTED_DATE* Fit_Laplacian_psi_1_3;
extern struct FUNC_FITTED_DATE* Fit_Laplacian_psi_1_4;


extern char* File_psi_1_0_out; // ψ_1(r)输出文件
extern char* File_psi_1_1_out;
extern char* File_psi_1_2_out;
extern char* File_psi_1_3_out;
extern char* File_psi_1_4_out;
extern char* File_psi_1_0_fit; // ψ_1(r)拟合结文件
extern char* File_psi_1_1_fit;
extern char* File_psi_1_2_fit;
extern char* File_psi_1_3_fit;
extern char* File_psi_1_4_fit;


extern char* File_Laplacian_psi_1_0_out; // ∆ψ_1(r)输出文件
extern char* File_Laplacian_psi_1_1_out;
extern char* File_Laplacian_psi_1_2_out;
extern char* File_Laplacian_psi_1_3_out;
extern char* File_Laplacian_psi_1_4_out;
extern char* File_Laplacian_psi_1_0_fit; // ∆ψ_1(r)拟合结文件
extern char* File_Laplacian_psi_1_1_fit;
extern char* File_Laplacian_psi_1_2_fit;
extern char* File_Laplacian_psi_1_3_fit;
extern char* File_Laplacian_psi_1_4_fit;


extern bool PT_profile_simplify; //是否启用简化版本的计算，typical profile 计算
extern bool PT_threshold_simplify;  //是否启用简化版本的计算，threshold 计算
extern bool PT_Mass_Relative; //计算黑洞的质量分布时，是否使用相对质量来进行表示和计算
extern bool PT_cal_r_m_fix; //在计算数密度时，是否重新计算 r_m
extern bool Transfer_Function; //是否加入转移函数
extern bool Continuum_spectrum_cal_simplify; //连续谱计算是否采用简化
extern bool Continuum_spectrum_judge_help; //连续谱功率谱判断辅助
extern bool Stdout_verbose; //命令行输出


//Press_Schechter 相关计算使用
extern arb_t PS_Int_variance_min; // 计算方差 XX YY XY 积分
extern arb_t PS_Int_variance_max;
extern arb_t PS_Int_variance_precision;

extern arb_t PS_abundance_beta_delta_k_all_min; //用δ谱求连续谱
extern arb_t PS_abundance_beta_delta_k_all_max;
extern arb_t PS_abundance_beta_delta_k_all_precision;
extern ulong PS_abundance_beta_delta_k_all_Int_iterate_min;
extern ulong PS_abundance_beta_delta_k_all_Int_iterate_max;

extern arb_t PS_abundance_beta_delta_k_M_precision;
extern ulong PS_abundance_beta_delta_k_M_Int_iterate_min;
extern ulong PS_abundance_beta_delta_k_M_Int_iterate_max;



extern arb_t PS_Sigma_XX; //方差 XX YY XY 的值
extern arb_t PS_Sigma_XY;
extern arb_t PS_Sigma_YX;
extern arb_t PS_Sigma_YY;

extern arb_t PS_Sigma_XX_save; //方差 XX YY XY 保存
extern arb_t PS_Sigma_XY_save;
extern arb_t PS_Sigma_YY_save;

extern arb_t PS_Int_P_C_l_min; // 计算 C_ℓ 的概率密度分布 P(C_l)
extern arb_t PS_Int_P_C_l_max;
extern arb_t PS_Int_P_C_l_precision;

extern arb_t PS_Root_C_l_to_Y_min; // PS中δ情况，计算P(C_ℓ)需反解Y
extern arb_t PS_Root_C_l_to_Y_max;
extern arb_t PS_Root_C_l_to_Y_precision;
extern slong PS_Root_C_l_to_Y_num;

extern arb_ptr P_normalization_coefficient; //概率归一化系数


extern arb_t PS_M_ratio_max; // I 型扰动质量比的最大值
extern arb_t PT_M_ratio_max;
extern arb_t PS_C_th; // 压缩函数的临界值 C_th
extern arb_t PS_C_l_th; // 线性压缩函数的临界值 C_th(r)
extern arb_t PS_abundance_f_all_precision; //最终占比f积分的精度


extern arb_t Continuum_spectrum_A_old; //连续谱计算相关
extern arb_t Continuum_spectrum_x_m; //求连续谱的特征模式用
extern arb_t Continuum_spectrum_k_ch; //连续谱的特征模式
extern enum GET_K_CH_TYPE Get_k_ch_type; //求连续谱特征模式方法


//other
extern char Path_save[PATH_MAX+1]; //获取存储相对路径
extern char Out_date_file[PATH_MAX+30]; //数据输出文件
extern char Out_fitted_file[PATH_MAX+30]; //拟合数据输出文件
extern char Out_picture_file[PATH_MAX+30]; //画图数据输出文件
extern int time_begin,time_end; //计时用


#endif // __PBHS_PHY_CONSTANT__  
