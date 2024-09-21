#include "02_power_spectra.h"
#include "../gw_func.h"
#include <stdlib.h> 

//方法二

//y积分下界
static int Int_GW_power_spectra_02_x_y_a(arb_t res, const arb_t x, void* param, const slong order, slong prec)
{
    arb_t t;
    arb_init(t);
    
    arb_sub_ui(t,x,1,prec); //y=|1-x|
    arb_abs(res,t); 
    
    arb_clear(t);
    return 0;
}
//y积分上界
static int Int_GW_power_spectra_02_x_y_b(arb_t res, const arb_t x, void* param, const slong order, slong prec)
{
    arb_add_ui(res,x,1,prec); //y=1+x
    return 0;
}

//积分函数
static int interior_GW_power_spectra_02(arb_t res, const arb_t v, const arb_t u, void* param, const slong order, slong prec)
{
    arb_t s,t,w,x,v_u,v_2,u_2,k_v,k_u;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(x);
    arb_init(v_u);
    arb_init(v_2);
    arb_init(u_2);
    arb_init(k_v);
    arb_init(k_u);
    
    //接收传入参数
    struct Func_transfer_parameter *func_eta_k; //这里不需要手动分配，只需将其指向传入的指针即可
    func_eta_k=param;
    
    arb_mul(x,func_eta_k->p_2,func_eta_k->p_1,prec); //x=k*η
    arb_mul(v_u,v,u,prec); //v*u
    arb_sqr(v_2,v,prec); //v^2
    arb_sqr(u_2,u,prec); //u^2
    arb_mul(k_v,func_eta_k->p_2,v,prec); //k*v
    arb_mul(k_u,func_eta_k->p_2,u,prec); //k*u
    
    //4 * [ (4*v^2-(1+v^2-u^2)^2) / (4*v*u) ]^2 * I^2(v,u,x) * P_ζ(kv) * P_ζ(ku)
    
    arb_sub(s,v_2,u_2,prec); //4*[ (4*v^2-(1+v^2-u^2)^2) / (4*v*u) ]^2
    arb_add_ui(s,s,1,prec);
    arb_sqr(s,s,prec);
    arb_mul_ui(t,v_2,4,prec);
    arb_sub(s,t,s,prec);
    
    arb_mul_ui(t,v_u,4,prec);
    arb_div(s,s,t,prec);
    arb_sqr(s,s,prec);
    arb_mul_ui(s,s,4,prec);
    
    I_func_v_u_x(t, v, u, x, prec); //I^2(v,u,x) * P_ζ(kv) * P_ζ(ku)
    arb_sqr(t,t,prec);
    arb_mul(s,s,t,prec);
    
    arb_log(w,k_v,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(s,s,t,prec);
    
    arb_log(w,k_u,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(x);
    arb_clear(v_u);
    arb_clear(v_2);
    arb_clear(u_2);
    arb_clear(k_v);
    arb_clear(k_u);
    
    return 0;
}

int GW_power_spectra_02(arb_t res, const arb_t eta, const arb_t k, slong prec)
{
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    //其中η和k是积分函数的参数，需传入
    //这里，对于结构体 Func_transfer_parameter 需手动分配内存
    struct Func_transfer_parameter *func_eta_k = (struct Func_transfer_parameter *)calloc(1,sizeof(struct Func_transfer_parameter));
    
    arb_init(func_eta_k->p_1); //使用arb_t变量前初始化
    arb_init(func_eta_k->p_2);
    
    arb_set(func_eta_k->p_1,eta); //设定传递参数
    arb_set(func_eta_k->p_2,k);
    
    //二元积分
    integration_binary_func(s, interior_GW_power_spectra_02, func_eta_k, 0,
                            Int_GW_power_spectra_min, Int_GW_power_spectra_max, Int_GW_power_spectra_precision, 
                            Int_GW_power_spectra_iterate_min, Int_GW_power_spectra_iterate_max,
                            Int_GW_power_spectra_02_x_y_a, NULL, 0,
                            Int_GW_power_spectra_02_x_y_b, NULL, 0,
                            Int_GW_power_spectra_precision,
                            Int_GW_power_spectra_iterate_min, Int_GW_power_spectra_iterate_max,
                            prec);
    
    arb_set(res,s);
    
    arb_clear(s);
    arb_clear(t);
    
    arb_clear(func_eta_k->p_1);
    arb_clear(func_eta_k->p_2);
    free(func_eta_k);
    
    return 0;
}


static int interior_GW_current_energy_density_02(arb_t res, const arb_t v, const arb_t u, void* param, const slong order, slong prec)
{
    arb_t s,t,w,q,v_u,v_2,u_2,k_v,k_u,sqrt_3;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(q);
    arb_init(v_u);
    arb_init(v_2);
    arb_init(u_2);
    arb_init(k_v);
    arb_init(k_u);
    arb_init(sqrt_3);
    
    //接收传入参数
    struct Func_transfer_parameter *func_k; //这里不需要手动分配，只需将其指向传入的指针即可
    func_k=param;
    
    arb_mul(v_u,v,u,prec); //v*u
    arb_sqr(v_2,v,prec); //v^2
    arb_sqr(u_2,u,prec); //u^2
    arb_mul(k_v,func_k->p_1,v,prec); //k*v
    arb_mul(k_u,func_k->p_1,u,prec); //k*u
    
    //具体表达式可见 2005.12306 (10)
    
    //3*T(v,u)
    arb_one(s); //3/4
    arb_mul_ui(s,s,3,prec);
    arb_sqrt(sqrt_3,s,prec); //√3
    arb_div_ui(s,s,4,prec);
    
    arb_sub(t,v_2,u_2,prec); //左边项
    arb_add_ui(t,t,1,prec);
    arb_sqr(t,t,prec);
    arb_mul_ui(w,v_2,4,prec);
    arb_sub(t,w,t,prec);
    arb_mul_ui(w,v_u,4,prec);
    arb_div(t,t,w,prec);
    arb_sqr(t,t,prec);
    arb_mul(s,s,t,prec);
    
    arb_add(t,v_2,u_2,prec); //中间项
    arb_sub_ui(t,t,3,prec);
    arb_mul_ui(w,v_u,2,prec);
    arb_div(t,t,w,prec);
    arb_pow_ui(t,t,4,prec);
    arb_mul(s,s,t,prec);
    
    
    arb_add(t,v,u,prec); //右边大括号 第一项，小括号第一项
    arb_sqr(t,t,prec);
    arb_neg(t,t);
    arb_add_ui(t,t,3,prec);
    arb_sub(w,u,v,prec);
    arb_sqr(w,w,prec);
    arb_neg(w,w);
    arb_add_ui(w,w,3,prec);
    arb_div(t,t,w,prec);
    arb_abs(t,t);
    arb_log(t,t,prec);
    
    arb_mul_ui(w,v_u,4,prec); //右边大括号 第一项，小括号第二项
    arb_add(q,u_2,v_2,prec);
    arb_sub_ui(q,q,3,prec);
    arb_div(w,w,q,prec);
    
    arb_sub(t,t,w,prec); //右边大括号，小括号内值
    arb_sqr(t,t,prec);
    
    arb_add(q,v,u,prec); //右边大括号 第二项
    arb_sub(q,q,sqrt_3,prec);
    Heaviside_Theta_function(w,q,prec);
    arb_const_pi(q,prec);
    arb_sqr(q,q,prec);
    arb_mul(w,w,q,prec);
    
    arb_add(t,t,w,prec); //右边大括号内两项相加
    
    arb_mul(s,s,t,prec); //3*T(v,u) 完
    
    arb_sqr(t,v_u,prec); //积分式左边
    arb_div(s,s,t,prec);
    
    arb_log(w,k_v,prec); //积分式右边功率谱
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(s,s,t,prec);
    
    arb_log(w,k_u,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(q);
    arb_clear(v_u);
    arb_clear(v_2);
    arb_clear(u_2);
    arb_clear(k_v);
    arb_clear(k_u);
    arb_clear(sqrt_3);
    
    return 0;
}

int GW_current_energy_density_02(arb_t res, const arb_t k, slong prec)
{
    arb_t s,t,w;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //其中η和k是积分函数的参数，需传入
    //这里，对于结构体 Func_transfer_parameter 需手动分配内存
    struct Func_transfer_parameter *func_k = (struct Func_transfer_parameter *)calloc(1,sizeof(struct Func_transfer_parameter));
    
    arb_init(func_k->p_1);//使用arb_t变量前初始化
    
    arb_set(func_k->p_1,k); //设定传递参数
    
    
    //积分前面的系数 c_g * Ω_{r,0}
    GW_spectra_convert_coefficient(s,prec);
    
    //二元积分
    integration_binary_func(w, interior_GW_current_energy_density_02, func_k, 0,
                            Int_GW_power_spectra_min, Int_GW_power_spectra_max, Int_GW_power_spectra_precision, 
                            Int_GW_power_spectra_iterate_min, Int_GW_power_spectra_iterate_max,
                            Int_GW_power_spectra_02_x_y_a, NULL, 0,
                            Int_GW_power_spectra_02_x_y_b, NULL, 0,
                            Int_GW_power_spectra_precision,
                            Int_GW_power_spectra_iterate_min, Int_GW_power_spectra_iterate_max,
                            prec);
    
    arb_mul(res,s,w,prec);//乘上前面的系数
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    
    arb_clear(func_k->p_1);
    free(func_k);
    
    return 0;
}
