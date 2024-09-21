#ifndef __PBHS_POWER_SPECTRA_NUMERIC_INFLATION_ODES_SLOVE__   /* Include guard */
#define __PBHS_POWER_SPECTRA_NUMERIC_INFLATION_ODES_SLOVE__

#include <arb.h> //高精度实数运算
#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"
#include "../other.h"

//此模块的内容，直接源自于 [mxhphy](https://github.com/mxhphy) 的 Python 程序

// dense output 的类型为 Inflation_dense_t, 初始化 Inflation_dense_init, 清理 Inflation_dense_clear
// dense output 插值拟合函数 Inflation_interp_fit_func_odes
// 扰动ODEs参数 Inflation_perturb_ODEs_param_t,
// 扰动ODEs参数: 初始化 Inflation_perturb_ODEs_param_init, 清理 Inflation_perturb_ODEs_param_clear
// 常微分方程求解 Inflation_ODEs_solver


typedef ODEs_DOP853_dense_t Inflation_dense_t; //定义一个新 dense output 类型，方便后面切换不同的ODEs求解方法

//下面对三个关于 dense output 的函数定义别名，方便后面切换不同的ODEs求解方法
//此间的 Inflation_use_dense_init_func 其参数和返回值需与要使用的 dense output 初始化一样
typedef ODEs_DOP853_dense_t (*Inflation_use_dense_init_func)(slong num, slong dim);
extern const Inflation_use_dense_init_func Inflation_dense_init;

typedef void (*Inflation_use_dense_clear_func)(ODEs_DOP853_dense_t dense_out, slong num, slong dim);
extern const Inflation_use_dense_clear_func Inflation_dense_clear;

typedef void (*Inflation_use_interp_fit_func_odes_func)(arb_t res, const arb_t x,
                                                        ODEs_DOP853_dense_t dense_out, const slong i_y, //i_y 表示微分方程解 y 中的第 i 个
                                                        slong prec);
extern const Inflation_use_interp_fit_func_odes_func Inflation_interp_fit_func_odes;

//再对常微分方程求解器使用名
typedef void (*Inflation_use_ODEs_solver_func)(arb_ptr y_end, my_odes_func func, const slong dim, void *param, const slong order, //常微分方程组函数
                                               const arb_t x_start, const arb_ptr y_start, //给定初始条件
                                               const arb_t x_end, //求出点 x_end 对应的函数值
                                               const arb_t error_abs, const arb_t error_rel, //绝对误差和相对误差 //误差为绝对精度
                                               const ODEs_DOP853_dense_t dense_out, //此参数可以为空
                                               slong prec);
extern const Inflation_use_ODEs_solver_func Inflation_ODEs_solver;



//此结构体用于暴胀中扰动微分方程的参数传递
//此方法避免使用全局变量，可适应于多线程
struct Inflation_perturb_ODEs_param_structure {
    Inflation_dense_t d_out; //参数一，背景方程解的dense output
    arb_t fourier_k; //参数二，求解所需的Fourier模式
};

typedef struct Inflation_perturb_ODEs_param_structure* Inflation_perturb_ODEs_param_t;


//暴胀中扰动微分方程参数传递的结构体，其初始化和清理函数
Inflation_perturb_ODEs_param_t Inflation_perturb_ODEs_param_init(Inflation_dense_t d_out, arb_t fourier_k);
void Inflation_perturb_ODEs_param_clear(Inflation_perturb_ODEs_param_t param_p);


// model parameters, 全局变量
extern arb_t Inf_Mpl,Inf_Mpl_2,Inf_Mpl_4,Inf_Mpl_6; //Planck质量极其幂次
extern arb_t Inf_Delta_V,Inf_V0,Inf_Phi_e;
extern arb_t Inf_Epsilon_1,Inf_sqrt_E1,Inf_Epsilon_2,Inf_sqrt_E2;
extern arb_t Inf_Eta_1,Inf_Eta_2;
extern arb_t Inf_Lambda,Inf_Delta_Phi_usr,Inf_Phi_s;

void Inflation_set_model_parameters(slong prec); //暴胀模型各种初始参数设定：势能参数、慢滚参数等

//微分方程计算所需函数
void Inflation_Smoothing_Step_Function(arb_t res, const arb_t x, const arb_t smoothness, slong prec);
void Inflation_S_prime(arb_t res, const arb_t x, const arb_t smoothness, slong prec);
void Inflation_S_pp(arb_t res, const arb_t x, const arb_t smoothness, slong prec);
void Inflation_V_phi(arb_t res, const arb_t phi, slong prec);
void Inflation_V_phi_p(arb_t res, const arb_t phi, slong prec);
void Inflation_V_phi_pp(arb_t res, const arb_t phi, slong prec);
void Inflation_m_eff(arb_t res, const arb_t phi_dot, const arb_t phi_ddot,
                     const arb_t H, const arb_t H_dot, const arb_t V_phi_pp, slong prec);


//背景微分方程
int Inflation_background_coupled_odes(arb_ptr yp, const arb_t x, const arb_ptr y, const slong dim,
                      void* param, const slong order, slong prec);

//扰动微分方程
int Inflation_perturbation_phi_odes(arb_ptr yp, const arb_t t, const arb_ptr y, const slong dim,
                               void* param, const slong order, slong prec);

//得到暴胀中，模式 k 进入视界的时间及此时的 N
void Inflation_get_time_k_enter(arb_t t, arb_t N, const arb_t k, const arb_t t_a, const arb_t t_b, //在区间 [t_a, t_b] 内找根
                                const arb_t a_i, const Inflation_dense_t d_out, // a_i 为尺度因子的初始值
                                slong prec);

//通过 t_ini, a_ini, fk 得到扰动方程的初始条件
void Inflation_get_perturb_oeds_initial_condition(arb_ptr y_start, const arb_t t_ini, const arb_t a_ini,
                                                  const arb_t fk, slong prec);

//利用在 x_end 处扰动方程的解 perturb_s, 及背景方程的 dense output d_out
//来求解 fk 模式的功率谱振幅
void Inflation_power_spectra_cal_at_fk(arb_t res, const arb_t fk, const arb_t x_end,
                                       const arb_ptr perturb_s, const Inflation_dense_t d_out,
                                       slong prec);

void Inflation_power_spectra_numeric_cal(slong prec); //功率谱求解

#endif // __PBHS_POWER_SPECTRA_NUMERIC_INFLATION_ODES_SLOVE__ 
