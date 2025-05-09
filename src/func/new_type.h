#ifndef __PBHS_NEW_TYPE__   /* Include guard */
#define __PBHS_NEW_TYPE__

#include <arb.h> //高精度实数运算
#include <sys/queue.h> // BSD 队列

//定义新的类型，用于各种情况的判断


//用来存储多项式系数拟合后得到的结果，相同区间，不用每次计算，加快速度
struct Interpolation_coe_structure {
    slong num;
    arb_poly_t* coe_poly;
};

typedef struct Interpolation_coe_structure* Interp_coe_t;


//用于微分方程的密度输出（dense output）
struct ODEs_DOPRI54_dense_output_structure {
    //用于四阶的多项式拟合
    //y(xn+θh)=r1+θ(r2+(1−θ)(r3+θ(r4+(1−θ)r5)))
    slong num_keep; //可容纳点的个数
    slong num_real; //真实存入点的个数
    slong dim; //微分方程组维数
    arb_ptr x;
    arb_ptr h;
    arb_ptr* r_1;
    arb_ptr* r_2;
    arb_ptr* r_3;
    arb_ptr* r_4;
    arb_ptr* r_5;
};
typedef struct ODEs_DOPRI54_dense_output_structure* ODEs_DOPRI54_dense_t;


struct ODEs_RFK45_dense_output_structure {
    //用于三阶的多项式拟合
    //y(xn+θh)=(1-θ)*y_i+θ*y_{i+1}+θ*(θ-1)*( (1-2*θ)(y_{i+1}-y_i)+(θ-1)*h*yp_i+θ*h*yp_{i+1} )
    slong num_keep; //可容纳点的个数
    slong num_real; //真实存入点的个数
    slong dim; //微分方程组维数
    arb_ptr x;
    arb_ptr h;
    arb_ptr* y_0;
    arb_ptr* y_1;
    arb_ptr* h_yp_0; //h*yp_0
    arb_ptr* h_yp_1; //h*yp_1
};
typedef struct ODEs_RFK45_dense_output_structure* ODEs_RFK45_dense_t;


struct ODEs_DOP853_dense_output_structure {
    //用于七阶的多项式拟合
    //参见 https://github.com/robclewley/pydstool/blob/master/PyDSTool/integrator/dop853.c
    slong num_keep; //可容纳点的个数
    slong num_real; //真实存入点的个数
    slong dim; //微分方程组维数
    arb_ptr x;
    arb_ptr h;
    arb_ptr* rcont1;
    arb_ptr* rcont2;
    arb_ptr* rcont3;
    arb_ptr* rcont4;
    arb_ptr* rcont5;
    arb_ptr* rcont6;
    arb_ptr* rcont7;
    arb_ptr* rcont8;
};
typedef struct ODEs_DOP853_dense_output_structure* ODEs_DOP853_dense_t;




//对多种窗口函数定义其类型
enum WINDOW_FUNC_TYPE
{
    Real_space_top_hat,
    Fourier_space_top_hat,
    Gaussian_hat,
    Natural_hat_xx,
    Natural_hat_yy
};


//找根方法定义
enum GET_K_CH_TYPE
{
    delta_approximation,
    horizon_re_enter
};


//找根方法定义
enum FIND_ROOT_ENUM_TYPE
{
    Brent_method,
    bisection_methode
};

//积分方法定义
enum INTEGRAL_ENUM_TYPE
{
    gauss_kronrod_iterate,
    gauss_kronrod_recursive,
    double_exponential
};

//诱导引力波方法定义
enum GW_INDUCED_ENUM_TYPE
{
    Li_gauss,
    Kohri_02,
    Espinosa_01
};

//β与f间的转换系数
enum BETA_TO_F_ENUM_TYPE
{
    beta_f_general_I,
    beta_f_general_II,
    beta_f_general_III
};

//用来进行计算的函数，如积分、找根等
typedef int (*my_calc_func)(arb_t out, const arb_t inp, void* param, const slong order, slong prec);

typedef int (*my_calc_func_binary)(arb_t out, const arb_t in_x, const arb_t in_y, void* param, const slong order, slong prec);

typedef int (*my_odes_func)(arb_ptr yp, const arb_t x, const arb_ptr y, const slong dim,
                             void* param, const slong order, slong prec);


typedef int (*func_NULL)(); //此类型后接的参数可任意


//功率谱类型定义
enum POWER_ENUM_TYPE
{
    delta_type,
    lognormal_type,
    power_law_type,
    box_type,
    broken_power_law_type,
    link_cmb_type,
    upward_step_spectra_type,
    numerical_cal_type
};

//求临界值的两种方法
enum PT_MU_TH_enum
{
    average_method_simple,
    average_method_new,
    q_parameter_method_simple,
    q_parameter_method_new
};


typedef enum POWER_ENUM_TYPE POWER; //功率谱


// ζ 类型定义
enum ZETA_ENUM_TYPE
{
    gaussian_type,
    power_expansion_type,
    exponential_tail_type,
    up_step_type,
    narrow_step_1_type,
    narrow_step_1_2_type
};

typedef enum ZETA_ENUM_TYPE ZETA; //功率谱


//求P(C_l)时，功率谱是δ类型时反解Y用
struct Find_root_delta_C_l_Y
{
    //传入求根，两个参数
    arb_t C_l, A;
};

//考虑所有k模式时，求β(k_M)用
struct Delta_Beta_M_All_k
{
    //传入k_M，一个参数
    arb_t k_M;
};

//求根时，作类型判断，对找C(r)的最大值特殊处理
enum Root_type
{
    Root_Normal,
    Root_C_Max
};

//函数传递多个参数用，例如积分时传多个传数进入
struct Func_transfer_parameter
{
    arb_t p_1;
    arb_t p_2;
    arb_t p_3;
};


//定义一个数据存储结构，存储读从文件中读取的数据
//可利用 BSD 提供的单链列表，其为一个动态数据结构，可以无限增长
struct FUNC_INPUT_POINT {
    arb_t x;
    arb_t y;
    SLIST_ENTRY(FUNC_INPUT_POINT) func_input_point_entries;
};

//#pragma pack(1) //内存里，各变量对齐的方式，默认4字节，修改为1字节
struct FUNC_FITTED_DATE {
    //拟合的过程是个三次函数 ax^3+bx^2+cx+d
    arb_t a;
    arb_t b;
    arb_t c;
    arb_t d;
    arb_t x; //相应区间的下限 [x_i, x_i+1] -> x_i ，用来判断对应的取值
    SLIST_ENTRY(FUNC_FITTED_DATE) func_fitted_date_entries;
};
//#pragma pack() // 释放内存对齐





#endif // __PBHS_NEW_TYPE__  
