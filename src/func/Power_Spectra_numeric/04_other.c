#include "Inflation_odes_slove.h"
#include <stdlib.h> 

//为密度输出及插值定义了别名，方便后面换用不同方法
const Inflation_use_dense_init_func Inflation_dense_init= ODEs_DOP853_dense_init;
const Inflation_use_dense_clear_func Inflation_dense_clear= ODEs_DOP853_dense_clear;
const Inflation_use_interp_fit_func_odes_func Inflation_interp_fit_func_odes= Interpolation_fit_func_odes_DOP853;

//为求解器定义别名
const Inflation_use_ODEs_solver_func Inflation_ODEs_solver=ODEs_DOP853;


//这里，定义了程序计算中所需的其它各种函数

//有效质量，求解 Sasaki-Mukhanov 方程时用
void Inflation_m_eff(arb_t res, const arb_t phi_dot, const arb_t phi_ddot,
                     const arb_t H, const arb_t H_dot, const arb_t V_phi_pp, slong prec)
{
    arb_t s,t,w;    
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //V_phi_pp - (1/Mpl)**2 * (3.*phi_dot**2. + (2*phi_dot*phi_ddot)/H - (H_dot*(phi_dot**2.))/(H**2))
    
    //(3.*phi_dot**2. + (2*phi_dot*phi_ddot)/H - (H_dot*(phi_dot**2.))/(H**2))
    arb_sqr(s,phi_dot,prec);
    arb_mul_ui(s,s,3,prec);
    
    arb_mul(t,phi_dot,phi_ddot,prec);
    arb_mul_si(t,t,2,prec);
    arb_div(t,t,H,prec);
    arb_add(s,s,t,prec);
    
    arb_sqr(t,phi_dot,prec);
    arb_mul(t,t,H_dot,prec);
    arb_sqr(w,H,prec);
    arb_div(t,t,w,prec);
    arb_sub(s,s,t,prec);
    
    arb_inv(t,Inf_Mpl,prec);
    arb_sqr(t,t,prec);
    arb_mul(s,s,t,prec);
    
    arb_sub(res,V_phi_pp,s,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
}


//暴胀中扰动微分方程参数传递的结构体，其初始化和清理函数
Inflation_perturb_ODEs_param_t Inflation_perturb_ODEs_param_init(Inflation_dense_t d_out, arb_t fourier_k)
{
    Inflation_perturb_ODEs_param_t param_p;
    param_p=(Inflation_perturb_ODEs_param_t)calloc(1,sizeof(struct Inflation_perturb_ODEs_param_structure));
    
    param_p->d_out=d_out;
    
    arb_init(param_p->fourier_k);
    arb_set(param_p->fourier_k,fourier_k);
    
    return param_p;
}

void Inflation_perturb_ODEs_param_clear(Inflation_perturb_ODEs_param_t param_p)
{
    arb_clear(param_p->fourier_k);
    free(param_p);
}


//得到暴胀中，模式 k 进入视界的时间及此时的 N
struct Get_time_k_enter_structure {
    Inflation_dense_t d_out;
    arb_t ln_a_i;
    arb_t ln_k;
};

static int interior_Inflation_get_time_k_enter(arb_t res, const arb_t ln_t,
                                   void * pass, const slong order,  slong prec)
{
    arb_t s,t,w,N,ln_H,phi,phi_dot,V;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(N);
    arb_init(ln_H);
    arb_init(phi);
    arb_init(phi_dot);
    arb_init(V);
    
    struct Get_time_k_enter_structure* param;
    param=pass;
    
    
    //进入视界条件 k=n*aH ⇒ ln(k)=N+ln(a_i)+ln(H)+ln(n)
    //利用这里得到的时间，作为后面计算的初始条件
    //系数 n 保证后面的计算是从视界内开始
    
    arb_exp(t,ln_t,prec); //恢复线性时间
    
    //背景解 N = N_interp(t)
    Inflation_interp_fit_func_odes(N, t, param->d_out, 3,  prec);
    
    //背景解 phi = phi_interp(t)
    Inflation_interp_fit_func_odes(phi, t, param->d_out, 1, prec);
    
    //背景解 phi_dot = phi_dot_interp(t)
    Inflation_interp_fit_func_odes(phi_dot, t, param->d_out, 0, prec);
    
    Inflation_V_phi(V,phi,0,prec); //原函数 V = V_phi(phi)
    
    //背景解 H = H_interp(t), 这里不用插值
    //H = np.sqrt((1./6.) * phi_dot**2 + V / 3.)
    arb_sqr(s,phi_dot,prec);
    arb_div_ui(s,s,6,prec);
    arb_div_ui(w,V,3,prec);
    arb_add(s,s,w,prec);
    arb_sqrt(ln_H,s,prec); //得到 H
    arb_log(ln_H,ln_H,prec); //取对数
    
    //ln(k)=N+ln(a_i)+ln(H)+ln(n) ⇒ N+ln(a_i)+ln(H)+ln(n)-ln(k)=0
    
    arb_add(s,N,param->ln_a_i,prec);
    arb_add(s,s,ln_H,prec);
    
    arb_set_si(w,100); // ln(n) 为了保证计算的准确性，一般取到 ln(50)
    arb_log(w,w,prec);
    arb_add(s,s,w,prec);
    
    arb_sub(res,s,param->ln_k,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(N);
    arb_clear(ln_H);
    arb_clear(phi);
    arb_clear(phi_dot);
    arb_clear(V);
    return 0;
}

//得到暴胀中，模式 k 进入视界的时间及此时的 N
//这里的 k 没取对数，以实际 k 输入
//这里，输出的时间 t 没取对数值
void Inflation_get_time_k_enter(arb_t t, arb_t N, const arb_t k, const arb_t t_a, const arb_t t_b, //在区间 [t_a, t_b] 内找根
                           const arb_t a_i, const Inflation_dense_t d_out, // a_i 为尺度因子的初始值
                           slong prec)
{
    //进入视界条件 k=aH ⇒ ln(k)=N+ln(a_i)+ln(H)
    arb_t s,ln_ta,ln_tb,error;
    
    arb_init(s);
    arb_init(ln_ta);
    arb_init(ln_tb);
    arb_init(error);
    
    
    struct Get_time_k_enter_structure* param; //用来传递参数
    
    param= (struct Get_time_k_enter_structure*)calloc(1,sizeof(struct Get_time_k_enter_structure)); //分配内存
    
    arb_init(param->ln_a_i); //使用前初始化其中变量
    arb_init(param->ln_k);
    param->d_out=d_out;
    
    arb_log(param->ln_a_i,a_i,prec);
    arb_log(param->ln_k,k,prec);
    
    //这里时间 t 会非常大，为了求根方便，我们用对数值计算
    if( arb_is_zero(t_a) ) //区间端点可能会有零的情况
    {
        arb_set_str(s,"1E-50",prec);
        arb_log(ln_ta,s,prec);
    }else
    {
        arb_log(ln_ta,t_a,prec);
    }
    
    if( arb_is_zero(t_b) )
    {
        arb_set_str(s,"1E-50",prec);
        arb_log(ln_tb,s,prec);
    }else
    {
        arb_log(ln_tb,t_b,prec);
    }
    
    
    
    arb_set_ui(error,1E5); //1E-5
    arb_inv(error,error,prec);
    
    //需要利用求根方法来得到 k 进入视界的时间
    Find_interval_root(s, interior_Inflation_get_time_k_enter, param, 0,
                       ln_ta, ln_tb, error,
                       10,Root_C_Max,
                       prec);
    arb_exp(t,s,prec);
    
    //求出此时的 N，可以通过此来计算 当前 a， ln(a)=N+ln(a_i) 
    Inflation_interp_fit_func_odes(s, t, param->d_out, 3,  prec);
    arb_set(N,s);
    
    arb_clear(s);
    arb_clear(ln_ta);
    arb_clear(ln_tb);
    arb_clear(error);
    
    arb_clear(param->ln_a_i); //释放结构体内存
    arb_clear(param->ln_k);
    free(param);
}


//通过 t_ini, a_ini, fk 得到扰动方程的初始条件
void Inflation_get_perturb_oeds_initial_condition(arb_ptr y_start, const arb_t t_ini, const arb_t a_ini,
                                             const arb_t fk, slong prec)
{
    arb_t s,t,cos_fka,sin_fka;
    arb_t Qphi_Real_t0, Qphi_Real_dot_t0, Qphi_Imag_t0, Qphi_Imag_dot_t0;
    arb_init(s);
    arb_init(t);
    arb_init(cos_fka);
    arb_init(sin_fka);
    
    arb_init(Qphi_Real_t0);
    arb_init(Qphi_Real_dot_t0);
    arb_init(Qphi_Imag_t0);
    arb_init(Qphi_Imag_dot_t0);
    
    
    //cos(f_k/a_i*t_i)
    arb_div(s,fk,a_ini,prec);
    arb_mul(s,s,t_ini,prec);
    arb_cos(cos_fka,s,prec);
    arb_sin(sin_fka,s,prec);//sin(f_k/a_i*t_i)
    
    // Initial conditions for perturbation
    //arb_one(a_ini); //a_ini = 1.0
    
    //Qphi_Real_t0 = cos(f_k/a_i*t_i) / (a_ini * np.sqrt(2. * fk))
    //Qphi_Imag_t0 = sin(f_k/a_i*t_i) / (a_ini * np.sqrt(2. * fk))
    arb_mul_si(s,fk,2,prec);
    arb_sqrt(s,s,prec);
    arb_mul(s,s,a_ini,prec);
    
    arb_div(Qphi_Real_t0,cos_fka,s,prec);
    arb_div(Qphi_Imag_t0,sin_fka,s,prec);
    
    //Qphi_Real_dot_t0 = -sin(f_k/a_i*t_i)*np.sqrt(fk) / (a_ini**2 * np.sqrt(2.))
    //Qphi_Imag_dot_t0 = cos(f_k/a_i*t_i)*np.sqrt(fk) / (a_ini**2 * np.sqrt(2.))
    arb_set_ui(s,2);
    arb_sqrt(s,s,prec);
    arb_sqr(t,a_ini,prec);
    arb_mul(t,t,s,prec);
    arb_sqrt(s,fk,prec);
    arb_div(s,s,t,prec);
    
    arb_mul(Qphi_Real_dot_t0,sin_fka,s,prec);
    arb_neg(Qphi_Real_dot_t0,Qphi_Real_dot_t0);
    arb_mul(Qphi_Imag_dot_t0,cos_fka,s,prec);
    
    // y_start = { Qphi_Real_t0, Qphi_Real_dot_t0, Qphi_Imag_t0, Qphi_Imag_dot_t0 }
    arb_set(y_start,Qphi_Real_t0);
    arb_set(y_start+1,Qphi_Real_dot_t0);
    arb_set(y_start+2,Qphi_Imag_t0);
    arb_set(y_start+3,Qphi_Imag_dot_t0);
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(cos_fka);
    arb_clear(sin_fka);
    
    arb_clear(Qphi_Real_t0);
    arb_clear(Qphi_Real_dot_t0);
    arb_clear(Qphi_Imag_t0);
    arb_clear(Qphi_Imag_dot_t0);
}



//利用在 x_end 处扰动方程的解 perturb_s, 及背景方程的 dense output d_out
//来求解 fk 模式的功率谱振幅
void Inflation_power_spectra_cal_at_fk(arb_t res, const arb_t fk, const arb_t x_end,
                            const arb_ptr perturb_s, const Inflation_dense_t d_out,
                            slong prec)
{
    arb_t s,t,w;
    arb_t Qphi_Real_end,Qphi_Imag_end;
    arb_t H_end,phi_dot_end,N_end,phi_end,V_end;
    
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(Qphi_Real_end);
    arb_init(Qphi_Imag_end);
    arb_init(H_end);
    arb_init(phi_dot_end);
    arb_init(N_end);
    arb_init(phi_end);
    arb_init(V_end);
    
    // Calculate P_R at the final time
    arb_set(Qphi_Real_end,perturb_s+0);
    arb_set(Qphi_Imag_end,perturb_s+2);
    
    //背景解 N = N_interp(t)
    Inflation_interp_fit_func_odes(N_end, x_end, d_out, 3,  prec);
    
    //背景解 phi = phi_interp(t)
    Inflation_interp_fit_func_odes(phi_end, x_end, d_out, 1, prec);
    
    //背景解 phi_dot = phi_dot_interp(t)
    Inflation_interp_fit_func_odes(phi_dot_end, x_end, d_out, 0, prec);
    
    Inflation_V_phi(V_end,phi_end,0,prec); //原函数 V = V_phi(phi)
    
    //背景解 H = H_interp(t), 这里不用插值
    //H = np.sqrt((1./6.) * phi_dot**2 + V / 3.)
    arb_sqr(s,phi_dot_end,prec);
    arb_div_ui(s,s,6,prec);
    arb_div_ui(w,V_end,3,prec);
    arb_add(s,s,w,prec);
    arb_sqrt(H_end,s,prec); //得到 H
    
    //P_R = ((fk**3) / (2 * np.pi**2)) * ((H_end*Qphi_Real_end/phi_dot_end)**2 + (H_end*Qphi_Imag_end/phi_dot_end)**2)
    
    //((fk**3) / (2 * np.pi**2))
    arb_sqr(t,Pi,prec);
    arb_mul_si(t,t,2,prec);
    arb_pow_ui(w,fk,3,prec);
    arb_div(w,w,t,prec);
    
    
    //(H_end*Qphi_Real_end/phi_dot_end)**2
    arb_mul(s,H_end,Qphi_Real_end,prec);
    arb_div(s,s,phi_dot_end,prec);
    arb_sqr(s,s,prec);
    
    //(H_end*Qphi_Imag_end/phi_dot_end)**2
    arb_mul(t,H_end,Qphi_Imag_end,prec);
    arb_div(t,t,phi_dot_end,prec);
    arb_sqr(t,t,prec);
    
    arb_add(s,s,t,prec);
    arb_mul(res,w,s,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(Qphi_Real_end);
    arb_clear(Qphi_Imag_end);
    arb_clear(H_end);
    arb_clear(phi_dot_end);
    arb_clear(N_end);
    arb_clear(phi_end);
    arb_clear(V_end);
}


//利用背景解，得到 H(t)
void Inflation_background_H_t(arb_t H, const arb_t t, const Inflation_dense_t d_out, slong prec)
{
    arb_t s,w,phi,phi_dot,V;
    arb_init(s);
    arb_init(w);
    arb_init(phi);
    arb_init(phi_dot);
    arb_init(V);
    
    Inflation_interp_fit_func_odes(phi, t, d_out, 1, prec); //ϕ
    Inflation_interp_fit_func_odes(phi_dot, t, d_out, 0, prec); //ϕ'
    Inflation_V_phi(V,phi,0,prec); //原函数 V = V_phi(phi)
    
    //背景解 H = H_interp(t), 这里不用插值
    //H = np.sqrt((1./6.) * phi_dot**2 + V / 3.)
    arb_sqr(s,phi_dot,prec);
    arb_div_ui(s,s,6,prec);
    arb_div_ui(w,V,3,prec);
    arb_add(s,s,w,prec);
    arb_sqrt(H,s,prec); //得到 H
    
    arb_clear(s);
    arb_clear(w);
    arb_clear(phi);
    arb_clear(phi_dot);
    arb_clear(V);
}


static int interior_phi_e_Inflation_get_model_g_h(arb_t res, const arb_t x,
                                               void * pass, const slong order,  slong prec)
{
    arb_t s;
    arb_init(s);
    
    Inflation_dense_t d_out; //密度输出通过参数传入
    d_out=pass;
    
    //求 ϕ_e 对应的时间
    Inflation_interp_fit_func_odes(s, x, d_out, 1, prec); //ϕ
    arb_sub(res,s,Inf_Phi_e,prec);
    
    arb_clear(s);
    return 0;
}

static int interior_phi_step_Inflation_get_model_g_h(arb_t res, const arb_t x,
                                                  void * pass, const slong order,  slong prec)
{
    arb_t s,w;
    arb_init(s);
    arb_init(w);
    
    Inflation_dense_t d_out; //密度输出通过参数传入
    d_out=pass;
    
    //求刚刚爬上 step  对应的时间
    Inflation_interp_fit_func_odes(s, x, d_out, 1, prec); //ϕ
    
    arb_inv(w,Inf_Lambda,prec); //刚爬上 step 的 ϕ = ϕ_e - 1/λ
    arb_sub(w,Inf_Phi_e,w,prec);
    
    arb_sub(res,s,w,prec);
    
    arb_clear(s);
    arb_clear(w);
    return 0;
}

struct interior_N_Inflation_get_model_g_h_s {
    arb_t N;
    Inflation_dense_t d_out;
};

static int interior_N_Inflation_get_model_g_h(arb_t res, const arb_t x,
                                                     void * pass, const slong order,  slong prec)
{
    arb_t s,w;
    arb_init(s);
    arb_init(w);
    
    struct interior_N_Inflation_get_model_g_h_s* p; //密度输出通过参数传入
    p=pass;
    
    //求刚刚爬上 step 后回归慢滚吸引子解对应的时间
    Inflation_interp_fit_func_odes(s, x, p->d_out, 3, prec); //N
    
    arb_sub(res,s,p->N,prec);
    
    arb_clear(s);
    arb_clear(w);
    return 0;
}

//输出相关信息，包括 h 和 g 等
void Inflation_get_model_correlated_info(const Inflation_dense_t d_out, slong prec)
{
    arb_t s,w,t_e,t_step,t_f,error,t_a,t_b,pi_c,pi_d,pi_f,h,g_cal,g_theory;
    arb_t N_e,N_step,epsilon_1,epsilon_2,eta_1,eta_2;
    arb_t nar_eta_1,nar_eta_2,nar_g,nar_kappa,nar_omega,nar_gamma,nar_beta,nar_A,nar_pi,nar_eq_h,amp_scale;
    
    arb_init(s);
    arb_init(w);
    arb_init(t_e);
    arb_init(t_step);
    arb_init(t_f);
    arb_init(error);
    arb_init(t_a);
    arb_init(t_b);
    arb_init(pi_c);
    arb_init(pi_d);
    arb_init(pi_f);
    arb_init(h);
    arb_init(g_cal);
    arb_init(g_theory);
    arb_init(N_e);
    arb_init(N_step);
    arb_init(epsilon_1);
    arb_init(epsilon_2);
    arb_init(eta_1);
    arb_init(eta_2);
    
    arb_init(nar_eta_1);
    arb_init(nar_eta_2);
    arb_init(nar_g);
    arb_init(nar_kappa);
    arb_init(nar_omega);
    arb_init(nar_gamma);
    arb_init(nar_beta);
    arb_init(nar_A);
    arb_init(nar_pi);
    arb_init(nar_eq_h);
    arb_init(amp_scale);
    
    arb_set_str(t_a,"1.6E6",prec);
    arb_set_str(t_b,"2E6",prec);
    
    arb_set_ui(error,1E14); //1E-5
    arb_inv(error,error,prec);
    
    //step 前的时间 t_e
    Find_interval_root(t_e, interior_phi_e_Inflation_get_model_g_h, d_out, 0,
                       t_a, t_b, error,
                       10,Root_Normal,
                       prec);
    //step 后的时间 t_step
    Find_interval_root(t_step, interior_phi_step_Inflation_get_model_g_h, d_out, 0,
                       t_a, t_b, error,
                       10,Root_Normal,
                       prec);
    
    //利用 N 估计回归慢滚吸引子时间
    Inflation_interp_fit_func_odes(N_e, t_e, d_out, 3, prec); //N
    Inflation_interp_fit_func_odes(N_step, t_step, d_out, 3, prec); //N
    arb_add_ui(w,N_step,2,prec); //step后差不多 1 e-fold 就回归到吸引子轨道，取 2 保险
    
    struct interior_N_Inflation_get_model_g_h_s* p;
    p=(struct interior_N_Inflation_get_model_g_h_s*)calloc(1,sizeof(struct interior_N_Inflation_get_model_g_h_s));
    
    arb_init(p->N);
    
    arb_set(p->N,w);
    p->d_out=d_out;
    
    //step 后回归慢滚吸引子解对应的时间 t_f (大概估计值)
    Find_interval_root(t_f, interior_N_Inflation_get_model_g_h, p, 0,
                       t_a, t_b, error,
                       10,Root_Normal,
                       prec);
    
    //step前的速度 π_c
    Inflation_interp_fit_func_odes(s, t_e, d_out, 0, prec); //ϕ'
    Inflation_background_H_t(w,t_e,d_out,prec); 
    arb_div(pi_c,s,w,prec); // dϕ/dN= ϕ'/H
    
    //step 后的速度 π_d
    Inflation_interp_fit_func_odes(s, t_step, d_out, 0, prec); //ϕ'
    Inflation_background_H_t(w,t_step,d_out,prec); 
    arb_div(pi_d,s,w,prec); // dϕ/dN= ϕ'/H
    
    //step 后回到慢滚吸引子解的速度 π_f
    Inflation_interp_fit_func_odes(s, t_f, d_out, 0, prec); //ϕ'
    Inflation_background_H_t(w,t_f,d_out,prec); 
    arb_div(pi_f,s,w,prec); // dϕ/dN= ϕ'/H
    
    //计算无量纲的慢滚参数
    //ϵ'_1=ϵ_1/V_0^2
    arb_sqr(w,Inf_V0,prec);
    arb_div(epsilon_1,Inf_Epsilon_1,w,prec);
    
    //η'_1=η_1/V_0^2, 注意这里与其它文章中的可能有个因子 2 的差异，根据定义表达式确定
    arb_sqr(w,Inf_V0,prec);
    arb_div(eta_1,Inf_Eta_1,w,prec);
    
    //ϵ'_2=ϵ_2/(V_0+ΔV)^2
    arb_add(w,Inf_V0,Inf_Delta_V,prec);
    arb_sqr(w,w,prec);
    arb_div(epsilon_2,Inf_Epsilon_2,w,prec);
    
    //η'_2=η_2/(V_0+ΔV), 注意这里与其它文章中的可能有个因子 2 的差异，根据定义表达式确定
    arb_add(w,Inf_V0,Inf_Delta_V,prec);
    arb_div(eta_2,Inf_Eta_2,w,prec);
    
    //数值计算所得 h 和 g
    arb_div(g_cal,pi_d,pi_c,prec); // g=π_d/π_c
    
    arb_mul_ui(w,epsilon_2,2,prec); // h≡-6*sqrt(2*ε_2)/π_d≈6 π_f/π_d, 其中的 ε_2 为无量纲的
    arb_sqrt(w,w,prec);
    arb_mul_ui(w,w,6,prec);
    arb_div(h,w,pi_d,prec);
    
    //理论计算所得 g
    arb_sqr(s,pi_c,prec); //π_d=sqrt(π_c^2-6*ΔV/V_0)
    arb_div(w,Inf_Delta_V,Inf_V0,prec);
    arb_mul_ui(w,w,6,prec);
    arb_sub(s,s,w,prec);
    arb_sqrt(s,s,prec);
    
    arb_div(g_theory,s,pi_c,prec); // g=π_d/π_c
    arb_abs(g_theory,g_theory);
    
    //
    //有限宽step对应非高斯性的相关参数
    //
    
    //η_i=2(2*ε_V_i-η_V_i)
    //其中的 ε_V 和 η_V 对应于上面的无量纲的 ε 和 η
    arb_mul_ui(s,epsilon_1,2,prec);
    arb_sub(s,s,eta_1,prec);
    arb_mul_ui(nar_eta_1,s,2,prec);
    
    arb_mul_ui(s,epsilon_2,2,prec);
    arb_sub(s,s,eta_2,prec);
    arb_mul_ui(nar_eta_2,s,2,prec);
    
    //g=π_2/π_2 与我们的定义一样
    arb_set(nar_g,g_cal);
    
    //κ=sqrt(ε_V_1/ε_V_2)
    arb_div(s,epsilon_1,epsilon_2,prec);
    arb_sqrt(nar_kappa,s,prec);
    
    //ω=ω_s2a≈sqrt(2)*|π_1|/Δϕ
    //Δϕ=1/λ
    arb_set_ui(s,2);
    arb_sqrt(s,s,prec);
    arb_abs(w,pi_c);
    arb_mul(s,s,w,prec);
    arb_mul(nar_omega,s,Inf_Lambda,prec);
    
    //利用 δN 形式计算非高斯性，算SR1的贡献时，需要个初始条件
    //由于 SR 的解为吸引子，我们设初始速度：π=π_c
    //这里，我们研究的模式主要是在step附近的，上面的做法问题不大
    //另外，其吸引子的速度约为 sqrt(2ε), 其为无量纲的慢滚参数，实际的值以数值计算结果为准
    arb_set(nar_pi,pi_c);
    
    //β=-1/pi
    arb_abs(s,nar_pi);
    arb_inv(nar_beta,s,prec);
    
    //γ= η_1*β/2 * (π/π_1)^(6/η_1) ≈ η_1*β/2 * (k/k_1)^3 , 其中 k 为 step 附近的模式
    arb_inv(s,nar_eta_1,prec);
    arb_mul_ui(s,s,6,prec);
    arb_div(w,nar_pi,pi_c,prec);
    arb_pow(s,w,s,prec);
    arb_mul(s,s,nar_beta,prec);
    arb_mul(s,s,nar_eta_1,prec);
    arb_div_ui(nar_gamma,s,2,prec);
    
    //A=β-κ*γ/(3*g)-γ/(ω*g^2)
    arb_mul_ui(w,nar_g,3,prec);
    arb_mul(s,nar_kappa,nar_gamma,prec);
    arb_div(s,s,w,prec);
    arb_sub(s,nar_beta,s,prec);
    arb_sqr(w,nar_g,prec);
    arb_mul(w,w,nar_omega,prec);
    arb_div(w,nar_gamma,w,prec);
    arb_sub(nar_A,s,w,prec);
    
    //对应的等效 h, h=2*γ/(A*g^2) 
    arb_sqr(s,nar_g,prec);
    arb_mul(s,s,nar_A,prec);
    arb_div(s,nar_gamma,s,prec);
    arb_mul_ui(nar_eq_h,s,2,prec);
    
    
    //估算放大的倍数 ε_1/ε_2/g^2
    arb_sqr(amp_scale,g_cal,prec);
    arb_mul(amp_scale,amp_scale,epsilon_2,prec);
    arb_div(amp_scale,epsilon_1,amp_scale,prec);
    
    
    printf("N_c: ");arb_printn(N_e, 20,0);printf("\n");
    printf("N_d: ");arb_printn(N_step, 20,0);printf("\n");
    printf("π_c: ");arb_printn(pi_c, 20,0);printf("\n");
    printf("π_d: ");arb_printn(pi_d, 20,0);printf("\n");
    printf("π_f: ");arb_printn(pi_f, 20,0);printf("\n\n");
    
    
    printf("h  : ");arb_printn(h, 20,0);printf("\n");
    printf("g_c: ");arb_printn(g_cal, 20,0);printf("\n\n");
    //printf("g_t: ");arb_printn(g_theory, 20,0);printf("\n\n");
    
    printf("ε_1: ");arb_printn(epsilon_1, 20,0);printf("\n");
    printf("ε_2: ");arb_printn(epsilon_2, 20,0);printf("\n");
    printf("η_1: ");arb_printn(eta_1, 20,0);printf("\n");
    printf("η_2: ");arb_printn(eta_2, 20,0);printf("\n\n");
    
    printf("narrow step β: ");arb_printn(nar_beta, 20,0);printf("\n");
    printf("narrow step κ: ");arb_printn(nar_kappa, 20,0);printf("\n");
    printf("narrow step g: ");arb_printn(nar_g, 20,0);printf("\n");
    printf("narrow step γ: ");arb_printn(nar_gamma, 20,0);printf("\n");
    printf("narrow step ω: ");arb_printn(nar_omega, 20,0);printf("\n");
    printf("narrow step A: ");arb_printn(nar_A, 20,0);printf("\n");
    printf("narrow step eq_h: ");arb_printn(nar_eq_h, 20,0);printf("\n\n");
    
    printf("amplify scale: ");arb_printn(amp_scale, 20,0);printf("\n\n");
    
    arb_clear(s);
    arb_clear(w);
    arb_clear(t_e);
    arb_clear(t_step);
    arb_clear(t_f);
    arb_clear(error);
    arb_clear(t_a);
    arb_clear(t_b);
    arb_clear(pi_c);
    arb_clear(pi_d);
    arb_clear(pi_f);
    arb_clear(h);
    arb_clear(g_cal);
    arb_clear(g_theory);
    arb_clear(N_e);
    arb_clear(N_step);
    arb_clear(epsilon_1);
    arb_clear(epsilon_2);
    arb_clear(eta_1);
    arb_clear(eta_2);
    
    arb_clear(nar_eta_1);
    arb_clear(nar_eta_2);
    arb_clear(nar_g);
    arb_clear(nar_kappa);
    arb_clear(nar_omega);
    arb_clear(nar_gamma);
    arb_clear(nar_beta);
    arb_clear(nar_A);
    arb_clear(nar_pi);
    arb_clear(nar_eq_h);
    arb_clear(amp_scale);
    
    arb_clear(p->N);
    free(p);
}
