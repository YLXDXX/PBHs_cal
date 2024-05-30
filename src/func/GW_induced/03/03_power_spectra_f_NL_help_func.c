#include "03_power_spectra.h"
#include <stdlib.h> 

static void y_ij(arb_t res, const arb_t t_i, const arb_t t_j,
                  const arb_t s_i, const arb_t s_j,
                  const arb_t phi_ij, slong prec)
{
    //具体表达式可见 2305.19950 (4.28)
    arb_t s,t,q;
    arb_init(s);
    arb_init(t);
    arb_init(q);
    
    //左边部分
    arb_add_ui(t,t_i,2,prec);
    arb_mul(t,t,t_i,prec);
    arb_sqr(s,s_i,prec);
    arb_neg(s,s);
    arb_add_ui(s,s,1,prec);
    arb_mul(t,t,s,prec);
    arb_mul(t,t,t_j,prec);
    
    arb_add_ui(s,t_j,2,prec);
    arb_mul(t,t,s,prec);
    arb_sqr(s,s_j,prec);
    arb_neg(s,s);
    arb_add_ui(s,s,1,prec);
    arb_mul(t,t,s,prec);
    arb_sqrt(t,t,prec);
    
    arb_cos(s,phi_ij,prec);
    arb_mul(t,t,s,prec);
    arb_div_ui(t,t,4,prec);
    
    //右边部分
    arb_add_ui(s,t_i,1,prec);
    arb_mul(s,s,s_i,prec);
    arb_neg(s,s);
    arb_add_ui(s,s,1,prec);
    arb_div_ui(s,s,4,prec);
    
    arb_add_ui(q,t_j,1,prec);
    arb_mul(q,q,s_j,prec);
    arb_neg(q,q);
    arb_add_ui(q,q,1,prec);
    arb_mul(s,s,q,prec);
    
    arb_add(res,t,s,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(q);
}

static void w_ij(arb_t res, const arb_t v_i, const arb_t v_j,
                  const arb_t t_i, const arb_t t_j,
                  const arb_t s_i, const arb_t s_j,
                  const arb_t phi_ij, slong prec)
{
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    //具体表达式可见 2305.19950 (4.29)
    arb_sqr(s,v_i,prec);
    arb_sqr(t,v_j,prec);
    arb_add(s,s,t,prec);
    
    y_ij(t, t_i, t_j, s_i, s_j, phi_ij, prec);
    
    arb_sub(s,s,t,prec);
    
    arb_sqrt(res,s,prec);
    
    arb_clear(s);
    arb_clear(t);
}

static void w_func_123(arb_t res, const arb_t v_1, const arb_t v_2, const arb_t v_3,
                  const arb_t t_1, const arb_t t_2, const arb_t t_3,
                  const arb_t s_1, const arb_t s_2, const arb_t s_3,
                  const arb_t phi_12, const arb_t phi_13, const arb_t phi_23, slong prec)
{
    //注意，这时 phi_ij 三个变量不独立 phi_12 + phi_23 = phi_13
    arb_t s,t,q;
    arb_init(s);
    arb_init(t);
    arb_init(q);
    
    //具体表达式可见 2305.19950 (4.30)
    arb_sqr(s,v_1,prec);
    arb_sqr(t,v_2,prec);
    arb_sqr(q,v_3,prec);
    arb_add(s,s,t,prec);
    arb_add(s,s,q,prec);
    
    y_ij(t, t_1, t_2, s_1, s_2, phi_12, prec); //y_12
    arb_add(s,s,t,prec);
    
    y_ij(t, t_1, t_3, s_1, s_3, phi_13, prec); //y_13
    arb_sub(s,s,t,prec);
    
    y_ij(t, t_2, t_3, s_2, s_3, phi_23, prec); //y_23
    arb_sub(s,s,t,prec);
    
    arb_sqrt(res,s,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(q);
}



int interior_GW_current_energy_density_Omega_H_4(arb_t res, const arb_t t_1, const arb_t s_1,
                                                      const arb_t t_2, const arb_t s_2,
                                                      void* param, const slong order, slong prec)
{
    arb_t s,t,w,u_1,v_1,u_2,v_2,k;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(u_1);
    arb_init(v_1);
    arb_init(u_2);
    arb_init(v_2);
    arb_init(k);
    
    //具体表达式可见 2305.19950 (4.21)
    
    //作了变量替换，变成矩形区域
    //u_1=(t_1 + s_1 + 1)/2
    arb_add(u_1,t_1,s_1,prec);//u_1
    arb_add_ui(u_1,u_1,1,prec);
    arb_div_ui(u_1,u_1,2,prec);
    
    arb_add(u_2,t_2,s_2,prec);//u_2
    arb_add_ui(u_2,u_2,1,prec);
    arb_div_ui(u_2,u_2,2,prec);
    
    //v_1=(t_1 - s_1 + 1)/2
    arb_sub(v_1,t_1,s_1,prec);//v_1
    arb_add_ui(v_1,v_1,1,prec);
    arb_div_ui(v_1,v_1,2,prec);
    
    arb_sub(v_2,t_2,s_2,prec);//v_2
    arb_add_ui(v_2,v_2,1,prec);
    arb_div_ui(v_2,v_2,2,prec);
    
    //接收传入参数
    struct Func_transfer_parameter *func_k; //这里不需要手动分配，只需将其指向传入的指针即可
    func_k=param;
    
    arb_set(k,func_k->p_1); //k
    
    
    //积分函数开始
    J_2_u_1_v_1_limited(s,u_1,v_1,prec);
    arb_mul(t,u_1,v_1,prec); 
    arb_mul(t,t,u_2,prec);
    arb_mul(t,t,v_2,prec);
    arb_sqr(t,t,prec);
    arb_div(s,s,t,prec);
    
    
    arb_mul(t,v_1,v_2,prec);
    arb_mul(t,t,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(s,s,t,prec);
    
    arb_mul(t,u_1,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(s,s,t,prec);
    
    
    arb_mul(t,v_1,u_2,prec);
    arb_mul(t,t,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(res,s,t,prec);
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(u_1);
    arb_clear(v_1);
    arb_clear(u_2);
    arb_clear(v_2);
    arb_clear(k);
    
    return 0;
}



int interior_GW_current_energy_density_Omega_R_6(arb_t res, const arb_t t_1, const arb_t s_1,
                                                      const arb_t t_2, const arb_t s_2,
                                                      const arb_t t_3, const arb_t s_3,
                                                      void* param, const slong order, slong prec)
{
    arb_t s,t,w,u_1,v_1,u_2,v_2,u_3,v_3,k;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(u_1);
    arb_init(v_1);
    arb_init(u_2);
    arb_init(v_2);
    arb_init(u_3);
    arb_init(v_3);
    arb_init(k);
    
    //具体表达式可见 2305.19950 (4.21)
    
    //作了变量替换，变成矩形区域
    //u_1=(t_1 + s_1 + 1)/2
    arb_add(u_1,t_1,s_1,prec);//u_1
    arb_add_ui(u_1,u_1,1,prec);
    arb_div_ui(u_1,u_1,2,prec);
    
    arb_add(u_2,t_2,s_2,prec);//u_2
    arb_add_ui(u_2,u_2,1,prec);
    arb_div_ui(u_2,u_2,2,prec);
    
    arb_add(u_3,t_3,s_3,prec);//u_3
    arb_add_ui(u_3,u_3,1,prec);
    arb_div_ui(u_3,u_3,2,prec);
    
    //v_1=(t_1 - s_1 + 1)/2
    arb_sub(v_1,t_1,s_1,prec);//v_1
    arb_add_ui(v_1,v_1,1,prec);
    arb_div_ui(v_1,v_1,2,prec);
    
    arb_sub(v_2,t_2,s_2,prec);//v_2
    arb_add_ui(v_2,v_2,1,prec);
    arb_div_ui(v_2,v_2,2,prec);
    
    arb_sub(v_3,t_3,s_3,prec);//v_3
    arb_add_ui(v_3,v_3,1,prec);
    arb_div_ui(v_3,v_3,2,prec);
    
    //接收传入参数
    struct Func_transfer_parameter *func_k; //这里不需要手动分配，只需将其指向传入的指针即可
    func_k=param;
    
    arb_set(k,func_k->p_1); //k
    
    
    //积分函数开始
    J_2_u_1_v_1_limited(s,u_1,v_1,prec);
    arb_mul(t,u_1,v_1,prec);
    arb_mul(t,t,u_2,prec);
    arb_mul(t,t,v_2,prec);
    arb_mul(t,t,u_3,prec);
    arb_mul(t,t,v_3,prec);
    arb_sqr(t,t,prec);
    arb_div(s,s,t,prec);
    
    
    arb_mul(t,v_1,v_2,prec);
    arb_mul(t,t,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(s,s,t,prec);
    
    arb_mul(t,v_1,u_2,prec);
    arb_mul(t,t,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(s,s,t,prec);
    
    arb_mul(t,u_1,v_3,prec);
    arb_mul(t,t,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(s,s,t,prec);
    
    
    arb_mul(t,u_1,u_3,prec);
    arb_mul(t,t,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(res,s,t,prec);
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(u_1);
    arb_clear(v_1);
    arb_clear(u_2);
    arb_clear(v_2);
    arb_clear(u_3);
    arb_clear(v_3);
    arb_clear(k);
    
    return 0;
}


int interior_GW_current_energy_density_Omega_C_5(arb_t res, const arb_t t_1, const arb_t s_1,
                                                      const arb_t t_2, const arb_t s_2, const arb_t phi_12,
                                                      void* param, const slong order, slong prec)
{
    arb_t s,t,w,u_1,v_1,u_2,v_2,k,w_12;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(u_1);
    arb_init(v_1);
    arb_init(u_2);
    arb_init(v_2);
    arb_init(k);
    arb_init(w_12);
    
    //具体表达式可见 2305.19950 (4.31)
    
    //作了变量替换，变成矩形区域
    //u_1=(t_1 + s_1 + 1)/2
    arb_add(u_1,t_1,s_1,prec);//u_1
    arb_add_ui(u_1,u_1,1,prec);
    arb_div_ui(u_1,u_1,2,prec);
    
    arb_add(u_2,t_2,s_2,prec);//u_2
    arb_add_ui(u_2,u_2,1,prec);
    arb_div_ui(u_2,u_2,2,prec);
    
    //v_1=(t_1 - s_1 + 1)/2
    arb_sub(v_1,t_1,s_1,prec);//v_1
    arb_add_ui(v_1,v_1,1,prec);
    arb_div_ui(v_1,v_1,2,prec);
    
    arb_sub(v_2,t_2,s_2,prec);//v_2
    arb_add_ui(v_2,v_2,1,prec);
    arb_div_ui(v_2,v_2,2,prec);
    
    //接收传入参数
    struct Func_transfer_parameter *func_k; //这里不需要手动分配，只需将其指向传入的指针即可
    func_k=param;
    
    arb_set(k,func_k->p_1); //k
    
    //积分函数开始
    arb_mul_ui(w,phi_12,2,prec);
    arb_cos(w,w,prec);
    J_u_1_v_1_and_u_2_v_2_limited(s,u_1,v_1,u_2,v_2,prec);
    arb_mul(s,s,w,prec);
    
    arb_mul(s,s,v_1,prec);//最前面的u_i * v_i, i ∈ {1,2}
    arb_mul(s,s,u_1,prec);
    arb_mul(s,s,v_2,prec);
    arb_mul(s,s,u_2,prec);
    
    //分母部分
    arb_pow_ui(t,v_2,3,prec);
    arb_pow_ui(w,u_2,3,prec);
    arb_mul(t,t,w,prec);
    w_ij(w_12, v_1, v_2, t_1, t_2, s_1, s_2, phi_12, prec); //w_12
    arb_pow_ui(w,w_12,3,prec);
    arb_mul(t,t,w,prec);
    
    arb_div(s,s,t,prec);
    
    //分子部分
    arb_mul(t,v_2,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(s,s,t,prec);
    
    arb_mul(t,u_2,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(s,s,t,prec);
    
    
    arb_mul(t,w_12,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(res,s,t,prec);
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(u_1);
    arb_clear(v_1);
    arb_clear(u_2);
    arb_clear(v_2);
    arb_clear(k);
    arb_clear(w_12);
    
    return 0;
}

int interior_GW_current_energy_density_Omega_Z_5(arb_t res, const arb_t t_1, const arb_t s_1,
                                                        const arb_t t_2, const arb_t s_2, const arb_t phi_12,
                                                        void* param, const slong order, slong prec)
{
    arb_t s,t,w,u_1,v_1,u_2,v_2,k,w_12;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(u_1);
    arb_init(v_1);
    arb_init(u_2);
    arb_init(v_2);
    arb_init(k);
    arb_init(w_12);
    
    //具体表达式可见 2305.19950 (4.31)
    
    //作了变量替换，变成矩形区域
    //u_1=(t_1 + s_1 + 1)/2
    arb_add(u_1,t_1,s_1,prec);//u_1
    arb_add_ui(u_1,u_1,1,prec);
    arb_div_ui(u_1,u_1,2,prec);
    
    arb_add(u_2,t_2,s_2,prec);//u_2
    arb_add_ui(u_2,u_2,1,prec);
    arb_div_ui(u_2,u_2,2,prec);
    
    //v_1=(t_1 - s_1 + 1)/2
    arb_sub(v_1,t_1,s_1,prec);//v_1
    arb_add_ui(v_1,v_1,1,prec);
    arb_div_ui(v_1,v_1,2,prec);
    
    arb_sub(v_2,t_2,s_2,prec);//v_2
    arb_add_ui(v_2,v_2,1,prec);
    arb_div_ui(v_2,v_2,2,prec);
    
    //接收传入参数
    struct Func_transfer_parameter *func_k; //这里不需要手动分配，只需将其指向传入的指针即可
    func_k=param;
    
    arb_set(k,func_k->p_1); //k
    
    //积分函数开始
    arb_mul_ui(w,phi_12,2,prec);
    arb_cos(w,w,prec);
    J_u_1_v_1_and_u_2_v_2_limited(s,u_1,v_1,u_2,v_2,prec);
    arb_mul(s,s,w,prec);
    
    arb_mul(s,s,v_1,prec);//最前面的u_i * v_i, i ∈ {1,2}
    arb_mul(s,s,u_1,prec);
    arb_mul(s,s,v_2,prec);
    arb_mul(s,s,u_2,prec);
     
    //分母部分
    arb_pow_ui(t,v_2,3,prec);
    arb_pow_ui(w,u_1,3,prec);
    arb_mul(t,t,w,prec);
    w_ij(w_12, v_1, v_2, t_1, t_2, s_1, s_2, phi_12, prec); //w_12
    arb_pow_ui(w,w_12,3,prec);
    arb_mul(t,t,w,prec);
    
    arb_div(s,s,t,prec);
    
    //分子部分
    arb_mul(t,v_2,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(s,s,t,prec);
    
    arb_mul(t,u_1,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(s,s,t,prec);
    
    
    arb_mul(t,w_12,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(res,s,t,prec);
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(u_1);
    arb_clear(v_1);
    arb_clear(u_2);
    arb_clear(v_2);
    arb_clear(k);
    arb_clear(w_12);
    
    return 0;
}

int interior_GW_current_energy_density_Omega_P_8(arb_t res, const arb_t t_1, const arb_t s_1,
                                                        const arb_t t_2, const arb_t s_2,
                                                        const arb_t t_3, const arb_t s_3,
                                                        const arb_t phi_12, const arb_t phi_23,
                                                        void* param, const slong order, slong prec)
{
    arb_t s,t,w,u_1,v_1,u_2,v_2,u_3,v_3,phi_13,w_13,w_23,k;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(u_1);
    arb_init(v_1);
    arb_init(u_2);
    arb_init(v_2);
    arb_init(u_3);
    arb_init(v_3);
    arb_init(phi_13);
    arb_init(w_13);
    arb_init(w_23);
    arb_init(k);
    
    //注意，这时 phi_ij 三个变量不独立 phi_12 + phi_23 = phi_13
    arb_add(phi_13,phi_12,phi_23,prec);
    
    
    //具体表达式可见 2305.19950 (4.33)
    
    //作了变量替换，变成矩形区域
    //u_1=(t_1 + s_1 + 1)/2
    arb_add(u_1,t_1,s_1,prec);//u_1
    arb_add_ui(u_1,u_1,1,prec);
    arb_div_ui(u_1,u_1,2,prec);
    
    arb_add(u_2,t_2,s_2,prec);//u_2
    arb_add_ui(u_2,u_2,1,prec);
    arb_div_ui(u_2,u_2,2,prec);
    
    arb_add(u_3,t_3,s_3,prec);//u_3
    arb_add_ui(u_3,u_3,1,prec);
    arb_div_ui(u_3,u_3,2,prec);
    
    //v_1=(t_1 - s_1 + 1)/2
    arb_sub(v_1,t_1,s_1,prec);//v_1
    arb_add_ui(v_1,v_1,1,prec);
    arb_div_ui(v_1,v_1,2,prec);
    
    arb_sub(v_2,t_2,s_2,prec);//v_2
    arb_add_ui(v_2,v_2,1,prec);
    arb_div_ui(v_2,v_2,2,prec);
    
    arb_sub(v_3,t_3,s_3,prec);//v_3
    arb_add_ui(v_3,v_3,1,prec);
    arb_div_ui(v_3,v_3,2,prec);
    
    //接收传入参数
    struct Func_transfer_parameter *func_k; //这里不需要手动分配，只需将其指向传入的指针即可
    func_k=param;
    
    arb_set(k,func_k->p_1); //k
    
    
    //积分函数开始
    arb_mul_ui(w,phi_12,2,prec);
    arb_cos(w,w,prec);
    J_u_1_v_1_and_u_2_v_2_limited(s,u_1,v_1,u_2,v_2,prec);
    arb_mul(s,s,w,prec);
    
    arb_mul(s,s,v_1,prec);//最前面的u_i * v_i, i ∈ {1,2,3}
    arb_mul(s,s,u_1,prec);
    arb_mul(s,s,v_2,prec);
    arb_mul(s,s,u_2,prec);
    arb_mul(s,s,v_3,prec);
    arb_mul(s,s,u_3,prec);
    
    //分母部分
    arb_pow_ui(t,v_3,3,prec);
    arb_pow_ui(w,u_3,3,prec);
    arb_mul(t,t,w,prec);
    
    w_ij(w_13, v_1, v_3, t_1, t_3, s_1, s_3, phi_13, prec); //w_13
    w_ij(w_23, v_2, v_3, t_2, t_3, s_2, s_3, phi_23, prec); //w_23
    
    arb_pow_ui(w,w_13,3,prec);
    arb_mul(t,t,w,prec);
    arb_pow_ui(w,w_23,3,prec);
    arb_mul(t,t,w,prec);
    
    arb_div(s,s,t,prec);
    
    //分子部分
    arb_mul(t,v_3,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(s,s,t,prec);
    
    arb_mul(t,u_3,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(s,s,t,prec);
    
    arb_mul(t,w_13,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(s,s,t,prec);
    
    arb_mul(t,w_23,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(res,s,t,prec);
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(u_1);
    arb_clear(v_1);
    arb_clear(u_2);
    arb_clear(v_2);
    arb_clear(u_3);
    arb_clear(v_3);
    arb_clear(phi_13);
    arb_clear(w_13);
    arb_clear(w_23);
    arb_clear(k);
    
    return 0;
}


int interior_GW_current_energy_density_Omega_N_8(arb_t res, const arb_t t_1, const arb_t s_1,
                                                        const arb_t t_2, const arb_t s_2,
                                                        const arb_t t_3, const arb_t s_3,
                                                        const arb_t phi_12, const arb_t phi_23,
                                                        void* param, const slong order, slong prec)
{
    arb_t s,t,w,u_1,v_1,u_2,v_2,u_3,v_3,phi_13,w_13,w_23,w_123,k;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(u_1);
    arb_init(v_1);
    arb_init(u_2);
    arb_init(v_2);
    arb_init(u_3);
    arb_init(v_3);
    arb_init(phi_13);
    arb_init(w_13);
    arb_init(w_23);
    arb_init(w_123);
    arb_init(k);
    
    //注意，这时 phi_ij 三个变量不独立 phi_12 + phi_23 = phi_13
    arb_add(phi_13,phi_12,phi_23,prec);
    
    
    //具体表达式可见 2305.19950 (4.34)
    
    //作了变量替换，变成矩形区域
    //u_1=(t_1 + s_1 + 1)/2
    arb_add(u_1,t_1,s_1,prec);//u_1
    arb_add_ui(u_1,u_1,1,prec);
    arb_div_ui(u_1,u_1,2,prec);
    
    arb_add(u_2,t_2,s_2,prec);//u_2
    arb_add_ui(u_2,u_2,1,prec);
    arb_div_ui(u_2,u_2,2,prec);
    
    arb_add(u_3,t_3,s_3,prec);//u_3
    arb_add_ui(u_3,u_3,1,prec);
    arb_div_ui(u_3,u_3,2,prec);
    
    //v_1=(t_1 - s_1 + 1)/2
    arb_sub(v_1,t_1,s_1,prec);//v_1
    arb_add_ui(v_1,v_1,1,prec);
    arb_div_ui(v_1,v_1,2,prec);
    
    arb_sub(v_2,t_2,s_2,prec);//v_2
    arb_add_ui(v_2,v_2,1,prec);
    arb_div_ui(v_2,v_2,2,prec);
    
    arb_sub(v_3,t_3,s_3,prec);//v_3
    arb_add_ui(v_3,v_3,1,prec);
    arb_div_ui(v_3,v_3,2,prec);
    
    //接收传入参数
    struct Func_transfer_parameter *func_k; //这里不需要手动分配，只需将其指向传入的指针即可
    func_k=param;
    
    arb_set(k,func_k->p_1); //k
    
    
    //积分函数开始
    arb_mul_ui(w,phi_12,2,prec);
    arb_cos(w,w,prec);
    J_u_1_v_1_and_u_2_v_2_limited(s,u_1,v_1,u_2,v_2,prec);
    arb_mul(s,s,w,prec);
    
    arb_mul(s,s,v_1,prec);//最前面的u_i * v_i, i ∈ {1,2,3}
    arb_mul(s,s,u_1,prec);
    arb_mul(s,s,v_2,prec);
    arb_mul(s,s,u_2,prec);
    arb_mul(s,s,v_3,prec);
    arb_mul(s,s,u_3,prec);
    
    //分母部分
    w_ij(w_13, v_1, v_3, t_1, t_3, s_1, s_3, phi_13, prec); //w_13
    w_ij(w_23, v_2, v_3, t_2, t_3, s_2, s_3, phi_23, prec); //w_23
    w_func_123(w_123, v_1, v_2, v_3, t_1, t_2, t_3, s_1, s_2, s_3, phi_12, phi_13, phi_23, prec); //w_123
    
    arb_pow_ui(t,u_3,3,prec);
    arb_pow_ui(w,w_13,3,prec);
    arb_mul(t,t,w,prec);
    arb_pow_ui(w,w_23,3,prec);
    arb_mul(t,t,w,prec);
    arb_pow_ui(w,w_123,3,prec);
    arb_mul(t,t,w,prec);
    
    arb_div(s,s,t,prec);
    
    //分子部分
    arb_mul(t,u_3,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(s,s,t,prec);
    
    arb_mul(t,w_13,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(s,s,t,prec);
    
    arb_mul(t,w_23,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(s,s,t,prec);
    
    arb_mul(t,w_123,k,prec);
    arb_log(w,t,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(res,s,t,prec);
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(u_1);
    arb_clear(v_1);
    arb_clear(u_2);
    arb_clear(v_2);
    arb_clear(u_3);
    arb_clear(v_3);
    arb_clear(phi_13);
    arb_clear(w_13);
    arb_clear(w_23);
    arb_clear(w_123);
    arb_clear(k);
    
    return 0;
}


