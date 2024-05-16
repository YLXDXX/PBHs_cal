#include "03_power_spectra.h"
#include <stdlib.h> 


static int interior_GW_current_energy_density_Omega_G(arb_t res, const arb_t t_1, const arb_t s_1, void* param, const slong order, slong prec)
{
    arb_t s,t,w,q,u_1,v_1,u_1_v_1,k_v_1,k_u_1;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(q);
    arb_init(u_1);
    arb_init(v_1);
    arb_init(u_1_v_1);
    arb_init(k_v_1);
    arb_init(k_u_1);
    
    //具体表达式可见 2305.19950 (4.21)
    
    //作了变量替换，变成矩形区域
    //u_1=(t_1 + s_1 + 1)/2
    arb_add(u_1,t_1,s_1,prec);
    arb_add_ui(u_1,u_1,1,prec);
    arb_div_ui(u_1,u_1,2,prec);
    
    //v_1=(t_1 - s_1 + 1)/2
    arb_sub(v_1,t_1,s_1,prec);
    arb_add_ui(v_1,v_1,1,prec);
    arb_div_ui(v_1,v_1,2,prec);
    
    //接收传入参数
    struct Func_transfer_parameter *func_k; //这里不需要手动分配，只需将其指向传入的指针即可
    func_k=param;
    
    arb_mul(u_1_v_1,u_1,v_1,prec); //u_1*v_1
    arb_mul(k_v_1,func_k->p_1,v_1,prec); //k*v_1
    arb_mul(k_u_1,func_k->p_1,u_1,prec); //k*u_1
    
   
    //积分函数开始
   J_2_u_1_v_1_limited(s,u_1,v_1,prec);
   arb_sqr(t,u_1_v_1,prec); 
   arb_div(s,s,t,prec);
   
    arb_log(w,k_u_1,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(s,s,t,prec);
    
    arb_log(w,k_v_1,prec);
    power_spectrum(t, w, prec); //这里传递的是以对数的形式传入
    arb_mul(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(q);
    arb_clear(u_1);
    arb_clear(v_1);
    arb_clear(u_1_v_1);
    arb_clear(k_v_1);
    arb_clear(k_u_1);
    
    return 0;
}

int GW_current_energy_density_Omega_G(arb_t res, const arb_t k, slong prec)
{
    arb_t s,t,w,c_g;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(c_g);
    
    //其中η和k是积分函数的参数，需传入
    //这里，对于结构体 Func_transfer_parameter 需手动分配内存
    struct Func_transfer_parameter *func_k = (struct Func_transfer_parameter *)calloc(1,sizeof(struct Func_transfer_parameter));
    
    arb_init(func_k->p_1);//使用arb_t变量前初始化
    
    arb_set(func_k->p_1,k); //设定传递参数
    
    
    //积分前面的系数 c_g * Ω_{r,0}
    arb_div(s,effective_g_star,effective_g_star_current,prec); //c_g= g_star/g_{star,0} * (g_{star,s,0}/g_{star,s})^{4/3}
    arb_div(t,effective_g_star_current_entropy,effective_g_star,prec);
    arb_one(w);
    arb_mul_ui(w,w,4,prec);
    arb_div_ui(w,w,3,prec);
    arb_pow(t,t,w,prec);
    arb_mul(c_g,s,t,prec);
    arb_mul(s,c_g,Omega_radiation,prec);
    
    arb_div_ui(s,s,3,prec); //(4.21)前积分系数
    
    //二元积分
    if( Int_GW_power_spectra_rectangle_adaptive ) //为真，采用自适应
    {
        integration_binary_rectangle_adaptive(w, interior_GW_current_energy_density_Omega_G, func_k, 0,
                                     Int_GW_power_spectra_x_min, Int_GW_power_spectra_x_max, Int_GW_power_spectra_x_precision, 
                                     Int_GW_power_spectra_iterate_x_min, Int_GW_power_spectra_iterate_x_max,
                                     Int_GW_power_spectra_y_min, Int_GW_power_spectra_y_max, Int_GW_power_spectra_y_precision, 
                                     Int_GW_power_spectra_iterate_y_min, Int_GW_power_spectra_iterate_y_max,
                                     prec);
    }else
    {
        integration_binary_rectangle(w, interior_GW_current_energy_density_Omega_G, func_k, 0,
                                     Int_GW_power_spectra_x_min, Int_GW_power_spectra_x_max, Int_GW_power_spectra_x_precision, 
                                     Int_GW_power_spectra_iterate_x_min, Int_GW_power_spectra_iterate_x_max,
                                     Int_GW_power_spectra_y_min, Int_GW_power_spectra_y_max, Int_GW_power_spectra_y_precision, 
                                     Int_GW_power_spectra_iterate_y_min, Int_GW_power_spectra_iterate_y_max,
                                     prec);
    }
    
    
    arb_mul(res,s,w,prec);//乘上前面的系数
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(c_g);
    free(func_k);
    
    return 0;
}

