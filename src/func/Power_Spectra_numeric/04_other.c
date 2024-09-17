#include "Inflation_odes_slove.h"
#include <stdlib.h> 

//为密度输出及插值定义了别名，方便后面换用不同方法
const Inflation_use_dense_init_func Inflation_dense_init= ODEs_DOPRI54_dense_init;
const Inflation_use_dense_clear_func Inflation_dense_clear= ODEs_DOPRI54_dense_clear;
const Inflation_use_interp_fit_func_odes_func Inflation_interp_fit_func_odes= Interpolation_fit_func_odes_DOPRI54;

//为求解器定义别名
const Inflation_use_ODEs_solver_func Inflation_ODEs_solver=ODEs_DOPRI54;


//这里，定义了程序计算中所需的其它函数

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
    
    
    //进入视界条件 k=aH ⇒ ln(k)=N+ln(a_i)+ln(H)
    arb_exp(t,ln_t,prec); //恢复线性时间
    
    //背景解 N = N_interp(t)
    Inflation_interp_fit_func_odes(N, t, param->d_out, 3,  prec);
    
    //背景解 phi = phi_interp(t)
    Inflation_interp_fit_func_odes(phi, t, param->d_out, 1, prec);
    
    //背景解 phi_dot = phi_dot_interp(t)
    Inflation_interp_fit_func_odes(phi_dot, t, param->d_out, 0, prec);
    
    Inflation_V_phi(V,phi,prec); //V = V_phi(phi)
    
    //背景解 H = H_interp(t), 这里不用插值
    //H = np.sqrt((1./6.) * phi_dot**2 + V / 3.)
    arb_sqr(s,phi_dot,prec);
    arb_div_ui(s,s,6,prec);
    arb_div_ui(w,V,3,prec);
    arb_add(s,s,w,prec);
    arb_sqrt(ln_H,s,prec); //得到 H
    arb_log(ln_H,ln_H,prec); //取对数
    
    //ln(k)=N+ln(a_i)+ln(H) ⇒ N+ln(a_i)+ln(H)-ln(k)=0
    
    arb_add(s,N,param->ln_a_i,prec);
    arb_add(s,s,ln_H,prec);
    
    //arb_set_str(w,"1.2",prec); //添加系数，将时间适当缩小
    //arb_mul(s,s,w,prec);
    
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
    
    Inflation_V_phi(V_end,phi_end,prec); //V = V_phi(phi)
    
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



