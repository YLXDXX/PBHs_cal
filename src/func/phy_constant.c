#include "phy_constant.h"
#include "new_type.h" //自定义的新类型
#include <arb.h> //高精度实数运算



//这里，用来定义各个全局变量、物理常数
//其它模快需要使用这里定义的变量
//只需要引用该模块就行 #include "phy_constant.h"


//一般变量，直接在此定义赋值即可
//double AA=3.141; 即可全局使用
//但是 arb_t 类型特殊
//在这里完成变量定义后，赋值操作需要在主程序中进行


//功率谱相关
POWER Power_spectrum_type; //功率谱类型
POWER Continuum_spectrum_type;
arb_t Power_A; //功率谱振幅
arb_t Power_sigma; //功率谱掌宽，例如 log-normal 的
arb_t Log_normal_mul_sigma; //功率谱掌宽的倍数
arb_t K_star; //功率谱参考尺度
arb_t Ln_K_star; //功率谱参考尺度，取对数
arb_t Box_K_middle; //Box功率谱 中心点
arb_t Box_width; //Box功率谱 宽度
arb_t BPL_alpha; //broken power law功率谱使用 α
arb_t BPL_beta; //broken power law功率谱使用 β
arb_t BPL_gamma; //broken power law功率谱使用 γ
arb_t Link_CMB_K_t; //link cmb 功率谱用
arb_t Link_CMB_K_m;
arb_t Link_CMB_A_s;
arb_t Link_CMB_A_b;
arb_t Link_CMB_K_star;
arb_t Link_CMB_P_m;
arb_t Link_CMB_n_s;
arb_t Link_CMB_n_b;
arb_t Upward_step_spectra_epsilon_S; //upward step模型的原生功率谱, SR-USR-SR
arb_t Upward_step_spectra_epsilon_V;
arb_t Upward_step_spectra_eta_V;
arb_t Upward_step_spectra_g;
arb_t Upward_step_spectra_h;
arb_t Upward_step_spectra_V_0;
arb_t Upward_step_spectra_Delta_V;
arb_t Upward_step_spectra_Delta_phi_USR;
arb_t Upward_step_spectra_phi_s;
arb_t Upward_step_spectra_phi_c;
arb_t Upward_step_spectra_k_s;
arb_t Upward_step_spectra_k_c;
arb_t Upward_step_spectra_ln_k_s;
arb_t Upward_step_spectra_ln_k_c;
arb_t Upward_step_spectra_tau_c;

//曲率扰动 ζ 相关
ZETA Zeta_type; // ζ 的类型
arb_t PT_mu; // ζ(r) 的参数
arb_t PT_mu_th; // ζ(r) 的参数 µ 的临界值
arb_t PT_mu_max; // ζ(r) 的参数 µ 的最大值
arb_t PT_k; // ζ(r) 的参数k
arb_t PT_beta_cal_need_exp_2_zeta_m; // 某个 M 对应的 e^(2*ζ_m)
arb_t R_MAX; // ζ(r) 取最大值时的 r 值
arb_t R_K_to_r_m; //求r_m动态区间用
arb_t Q_parameter_th; //ζ(r)取临界值时的q参数
arb_t R_m_times_K; //对于δ谱，r_m*k_*为一定值


//计算ζ(r)的辅助变量
arb_t Sigma_0_square; //(σ_n)^2
arb_t Sigma_1_square;
arb_t Sigma_2_square;
arb_t Sigma_3_square;
arb_t Sigma_4_square;
arb_t Gamma_1; //(γ_1)
arb_t Gamma_3; //(γ_3)
arb_t R_1; //(R_1)
arb_t R_3; //(R_3)

bool Peak_theory_source_zeta_gradient; //peak theory计算中，统计量是否取ζ的梯度




// C(r) 相关
ulong C_m_average_iterate_min; //求 C_m_average
ulong C_m_average_iterate_max;
arb_t C_m_average_precision; 


enum PT_MU_TH_enum PT_Mu_th_METHOD; //求临界值的方法
enum WINDOW_FUNC_TYPE PS_simple_window_func; //简单PS方法中窗口函数类型


//在视界进入时，视界质量M_H，形成的黑洞的质量为M，两者间的关系
//与曲率扰动 K 有关，与 γ 也有关
arb_t Mass_K;
arb_t Mass_gamma;
bool Critical_Collapse_Effect; //是否考虑临界坍缩效应


//宇宙学基本常数
arb_t Omega_DM; //暗物质所占比例
arb_t Omega_M; //baryons + dark matter 所占比例
arb_t Omega_B; //普通重子物质 baryons
arb_t Omega_Lambda; //暗能量 Ω_Λ
arb_t Omega_radiation; //辐射所占的比例 radiation density
arb_t Omega_gamma; //光子所占的比例 photon density
arb_t Hubble_constant_h; //哈勃常数 h
arb_t Hubble_constant_H_0; //哈勃常数 H_0
arb_t Scale_factor_a; // 尺度因子 a(t)，某些公式的需要
arb_t Scale_factor_a_0; //尺度因子 a(t_0)
arb_t Equation_Of_State_w; //状态方程 w
arb_t K_scale_eq; //在物质辐射相等时刻，对应的共动波数
arb_t Z_scale_eq; //在物质辐射相等时刻，对应的红移
arb_t A_scale_eq; //在物质辐射相等时刻，对应的尺度因子
arb_t T_scale_eq; //在物质辐射相等时刻，对应的温度

arb_t effective_g_star_eq; //在物质辐射相等时刻，对应的相对论有效自由度数目
arb_t effective_g_star_eq_entropy; //在物质辐射相等时刻，对应的熵有效自由度数目
arb_t effective_g_star; //重新进入视界后，对应的相对论有效自由度数目
arb_t effective_g_star_entropy; //重新进入视界后，对应的熵有效自由度数目
arb_t effective_g_star_current; //当今对应的相对论有效自由度数目
arb_t effective_g_star_current_entropy; //当今对应的熵有效自由度数目


//设置求值区间，如求 r_m PT_mu_th
arb_t Int_r_min;
arb_t Int_r_max;
slong Root_r_num;
arb_t Int_r_precision;
arb_t  Int_mu_min;
arb_t  Int_mu_max;
slong  Root_mu_num;
arb_t Int_mu_precision;


arb_t Int_sigma_n_min; //sigma_n
arb_t Int_sigma_n_max;
arb_t Int_sigma_n_precision;

arb_t Int_n_pk_k_min; // n_pk(mu,k) 中 k 的积分区间
arb_t Int_n_pk_k_max;
arb_t Int_n_pk_k_precision;



arb_t Root_M_to_mu_min; // M -> μ 求根用
arb_t Root_M_to_mu_max;
slong Root_M_to_mu_num;
arb_t Root_M_to_mu_precision;


ulong Integration_iterate_min; //积分的最小迭代次数
ulong Integration_iterate_max; //积分的最大迭代次数

// β 与 f 间的转换系数的选择
enum BETA_TO_F_ENUM_TYPE beta_to_f_type;


//非高斯相关参数
arb_t Power_expansion_f; //power-series expansion 二次项 f_NL -> A
arb_t Power_expansion_g; //power-series expansion 三次项 g_NL -> B
arb_t Power_expansion_four; //power-series expansion 四次项 four -> C
arb_t Power_expansion_five; //power-series expansion 五次项 five -> D
arb_t Power_expansion_six; //power-series expansion 六次项 six -> E
arb_t Up_step_h; //up-step 模型中用
arb_t Exponential_tail_beta; //exponential_tail 模型中用
arb_t Narrow_up_step_beta; // 有限宽 upward step 模型用
arb_t Narrow_up_step_kappa;
arb_t Narrow_up_step_g;
arb_t Narrow_up_step_gamma;
arb_t Narrow_up_step_omega;
arb_t Narrow_up_step_A;
arb_t Narrow_up_step_cutoff_1;

//诱导引力波
arb_t Int_GW_I_func_min; // I/I_c/I_s 积分用
arb_t Int_GW_I_func_max;
arb_t Int_GW_I_func_precision;
ulong Int_GW_I_func_iterate_min;
ulong Int_GW_I_func_iterate_max;

arb_t Int_GW_power_spectra_min; // 功率谱积分用
arb_t Int_GW_power_spectra_max;
arb_t Int_GW_power_spectra_precision;
ulong Int_GW_power_spectra_iterate_min;
ulong Int_GW_power_spectra_iterate_max;


arb_t Int_GW_power_spectra_x_min; // 功率谱积矩形分用
arb_t Int_GW_power_spectra_x_max;
arb_t Int_GW_power_spectra_x_precision;
ulong Int_GW_power_spectra_iterate_x_min;
ulong Int_GW_power_spectra_iterate_x_max;
arb_t Int_GW_power_spectra_y_min;
arb_t Int_GW_power_spectra_y_max;
arb_t Int_GW_power_spectra_y_precision;
ulong Int_GW_power_spectra_iterate_y_min;
ulong Int_GW_power_spectra_iterate_y_max;

bool Int_GW_power_spectra_rectangle_adaptive; //矩形二维积分是否采用自适应
enum GW_INDUCED_ENUM_TYPE GW_induced_method; //诱导引力波计算方法

bool GW_dim_integral_res_print;//积分结果打印
long int GW_dim_8_MINEVAL; //诱导引力波非高斯性积分参数
long int GW_dim_8_MAXEVAL;
long int GW_dim_8_NSTART;
long int GW_dim_8_NINCREASE;
long double GW_dim_8_EPSREL;
long double GW_dim_8_EPSABS;
char* GW_dim_8_t_upper;

long int GW_dim_6_MINEVAL;
long int GW_dim_6_MAXEVAL;
long int GW_dim_6_NSTART;
long int GW_dim_6_NINCREASE;
long double GW_dim_6_EPSREL;
long double GW_dim_6_EPSABS;
char* GW_dim_6_t_upper;

long int GW_dim_5_MINEVAL;
long int GW_dim_5_MAXEVAL;
long int GW_dim_5_NSTART;
long int GW_dim_5_NINCREASE;
long double GW_dim_5_EPSREL;
long double GW_dim_5_EPSABS;
char* GW_dim_5_t_upper;

long int GW_dim_4_MINEVAL;
long int GW_dim_4_MAXEVAL;
long int GW_dim_4_NSTART;
long int GW_dim_4_NINCREASE;
long double GW_dim_4_EPSREL;
long double GW_dim_4_EPSABS;
char* GW_dim_4_t_upper;

long int GW_dim_2_MINEVAL;
long int GW_dim_2_MAXEVAL;
long int GW_dim_2_NSTART;
long int GW_dim_2_NINCREASE;
long double GW_dim_2_EPSREL;
long double GW_dim_2_EPSABS;
char* GW_dim_2_t_upper;

//数学计算
arb_t Pi; //常数π
arb_t Pi_2; //常数2π (πx2)
slong Multithreaded_divide_integration_interval_number; //多线程计算区间分隔数目
slong Multithreaded_number; //多线程计算线程数

//slong prec; //控制计算精度
//slong N_n; //整数n，用于计算辅助

//arb_t _TEMP_arb_;//全局临时传递使用 arb_t 
//slong _TEMP_slong_;//全局临时传递使用 slong



//函数数据的输入输出，及拟合

arb_t Func_output_x_min; //输出范围
arb_t Func_output_x_max; 
ulong Func_output_number; //输出点的数目
char* Func_output_file ; //输出文件名
char* Func_output_delimiter; //输出数据间的分格符
unsigned int Func_output_skip_header; //输出中添加头开始几行的信息说明
unsigned int Func_output_column_x; //可能有多列数据，表明位置
unsigned int Func_output_column_y;
char * Func_output_fitted_file; //利用输出数据，拟合结果所存储名
size_t Fit_entries_number; //拟合读取数据的点数

//函数拟合相关
bool FIT_FUNC_IF; //是否开启拟合
struct FUNC_FITTED_DATE* Fit_test; //测试用

struct FUNC_FITTED_DATE* Fit_psi_1_0; // ψ_1(r) 及其各阶导数 拟合结果
struct FUNC_FITTED_DATE* Fit_psi_1_1;
struct FUNC_FITTED_DATE* Fit_psi_1_2;
struct FUNC_FITTED_DATE* Fit_psi_1_3;
struct FUNC_FITTED_DATE* Fit_psi_1_4;

struct FUNC_FITTED_DATE* Fit_Laplacian_psi_1_0; // ∆ψ_1(r) 及其各阶导数 拟合结果
struct FUNC_FITTED_DATE* Fit_Laplacian_psi_1_1;
struct FUNC_FITTED_DATE* Fit_Laplacian_psi_1_2;
struct FUNC_FITTED_DATE* Fit_Laplacian_psi_1_3;
struct FUNC_FITTED_DATE* Fit_Laplacian_psi_1_4;


char* File_psi_1_0_out; // ψ_1(r)输出文件
char* File_psi_1_1_out;
char* File_psi_1_2_out;
char* File_psi_1_3_out;
char* File_psi_1_4_out;
char* File_psi_1_0_fit; // ψ_1(r)拟合结文件
char* File_psi_1_1_fit;
char* File_psi_1_2_fit;
char* File_psi_1_3_fit;
char* File_psi_1_4_fit;


char* File_Laplacian_psi_1_0_out; // ∆ψ_1(r)输出文件
char* File_Laplacian_psi_1_1_out;
char* File_Laplacian_psi_1_2_out;
char* File_Laplacian_psi_1_3_out;
char* File_Laplacian_psi_1_4_out;
char* File_Laplacian_psi_1_0_fit; // ∆ψ_1(r)拟合结文件
char* File_Laplacian_psi_1_1_fit;
char* File_Laplacian_psi_1_2_fit;
char* File_Laplacian_psi_1_3_fit;
char* File_Laplacian_psi_1_4_fit;


bool PT_profile_simplify; //是否启用简化版本的计算，typical profile 计算
bool PT_threshold_simplify;  //是否启用简化版本的计算，threshold 计算
bool PT_Mass_Relative; //计算黑洞的质量分布时，是否使用相对质量来进行表示和计算
bool PT_cal_r_m_fix; //在计算数密度时，是否重新计算 r_m
bool Transfer_Function; //是否加入转移函数
bool Continuum_spectrum_cal_simplify; //连续谱计算是否采用简化
bool Continuum_spectrum_judge_help; //连续谱功率谱判断辅助
bool Stdout_verbose; //命令行输出

//Press_Schechter 相关计算使用
arb_t PS_Int_variance_min; // 计算方差 XX YY XY 积分
arb_t PS_Int_variance_max;
arb_t PS_Int_variance_precision;

arb_t PS_abundance_beta_delta_k_all_min; //用δ谱求连续谱
arb_t PS_abundance_beta_delta_k_all_max;
arb_t PS_abundance_beta_delta_k_all_precision;
ulong PS_abundance_beta_delta_k_all_Int_iterate_min;
ulong PS_abundance_beta_delta_k_all_Int_iterate_max;

arb_t PS_abundance_beta_delta_k_M_precision;
ulong PS_abundance_beta_delta_k_M_Int_iterate_min;
ulong PS_abundance_beta_delta_k_M_Int_iterate_max;



arb_t PS_Sigma_XX; //方差 XX YY XY 的值
arb_t PS_Sigma_XY;
arb_t PS_Sigma_YX;
arb_t PS_Sigma_YY;

arb_t PS_Sigma_XX_save; //方差 XX YY XY 保存
arb_t PS_Sigma_XY_save;
arb_t PS_Sigma_YY_save;

arb_t PS_Int_P_C_l_min; // 计算 C_ℓ 的概率密度分布 P(C_l)
arb_t PS_Int_P_C_l_max;
arb_t PS_Int_P_C_l_precision;

arb_t PS_Root_C_l_to_Y_min; // PS中δ情况，计算P(C_ℓ)需反解Y
arb_t PS_Root_C_l_to_Y_max;
arb_t PS_Root_C_l_to_Y_precision;
slong PS_Root_C_l_to_Y_num;

arb_t PS_Root_zeta_to_zeta_G_min;//计算概率 P(ζ), 需要反解 ζ= F(ζ_G)，有时反函数的解析表达式不能写出，且反函数不是单值
arb_t PS_Root_zeta_to_zeta_G_max;
arb_t PS_Root_zeta_to_zeta_G_precision;
slong PS_Root_zeta_to_zeta_G_num;

arb_ptr P_normalization_coefficient; //概率归一化系数

arb_t PS_M_ratio_max; // I 型扰动质量比的最大值
arb_t PT_M_ratio_max;
arb_t PS_M_ratio_min; // I 型扰动质量比的最小值
arb_t PT_M_ratio_min;
arb_t PS_C_th; // 压缩函数的临界值 C_th(r)
arb_t PS_C_l_th; // 线性压缩函数的临界值 C_th(r)
arb_t PS_zeta_th; // 曲率扰动估算阈值 ζ_th
arb_t PS_delta_th; // 密度扰动估算阈值 δ_th
arb_t PS_abundance_int_precision; //计算丰度β和f的积分精度
arb_t PS_abundance_simple_int_precision;
arb_t PS_abundance_simple_int_min;
arb_t PS_abundance_simple_int_max;


arb_t Continuum_spectrum_A_old; //连续谱计算相关
arb_t Delta_spectrum_x_m; //求连续谱的特征模式用
arb_t Delta_spectrum_reenter_coefficient_C; //求连续谱的特征模式用
arb_t Continuum_spectrum_k_ch; //连续谱的特征模式
arb_t Continuum_spectrum_k_ch_times_r_m; //k_ch*r_m
enum GET_K_CH_TYPE Get_k_ch_type; //求连续谱特征模式方法



//非高斯情况下，利用非高斯功率谱对 typical profile 修正用
bool Non_Gaussian_typical_profile_correction;


//other
char Path_save[PATH_MAX+1]=""; //获取存储相对路径
char Out_date_file[PATH_MAX+30]=""; //数据输出文件
char Out_fitted_file[PATH_MAX+30]=""; //拟合数据输出文件
char Out_fitted_x[PATH_MAX+30]=""; //用于拟合数函数
char Out_fitted_y[PATH_MAX+30]="";
char Out_fitted_NG_f_nl_k[PATH_MAX+30]=""; //用于拟合函数，非高斯性功率谱 f_NL 修正
char Out_fitted_NG_f_nl_P[PATH_MAX+30]="";
char Out_picture_file[PATH_MAX+30]=""; //画图数据输出文件

int time_begin,time_end; //计时用

arb_ptr FITTED_x; //用于拟合数函数
arb_ptr FITTED_y;
slong FITTED_num;
Interp_coe_t FITTED_interp_coe;

arb_ptr FITTED_NG_f_nl_k;  //用于拟合函数，非高斯性功率谱 f_NL 修正
arb_ptr FITTED_NG_f_nl_P;
slong FITTED_NG_f_nl_num;
Interp_coe_t FITTED_NG_f_nl_interp_coe;

