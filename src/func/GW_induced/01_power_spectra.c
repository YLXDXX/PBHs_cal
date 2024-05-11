#include "01_power_spectra.h"
#include <stdlib.h> 

//功率谱积分下界函数
static int Int_GW_power_spectra_01_x_y_a(arb_t res, const arb_t x, void* param, const slong order, slong prec)
{
    arb_t t;
    arb_init(t);
    
    arb_one(t);
    if( arb_lt(x,t) ) //x<1
    {
        arb_neg(t,x);//y=-x+1
        arb_add_ui(res,t,1,prec);
    }else //x≥1
    {
        arb_sub_ui(res,x,1,prec);//y=x-1
    }
    
    arb_clear(t);
    return 0;
}


//功率谱积分上界函数
static int Int_GW_power_spectra_01_x_y_b(arb_t res, const arb_t x, void* param, const slong order, slong prec)
{
    arb_add_ui(res,x,1,prec); //y=x+1
    return 0;
}

//功率谱积分函数f(x,y)
static int interior_GW_power_spectra_01_x_y(arb_t res, const arb_t x, const arb_t y, void* param, const slong order, slong prec)
{
    arb_t s,t,w,q,k_x,k_y,k_eta,x_2,y_2,I_c,I_s;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(q);
    arb_init(k_x);
    arb_init(k_y);
    arb_init(k_eta);
    arb_init(x_2);
    arb_init(y_2);
    arb_init(I_c);
    arb_init(I_s);
    
    //接收传入参数
    struct Func_transfer_parameter *func_eta_k; //这里不需要手动分配，只需将其指向传入的指针即可
    func_eta_k=param;
    
    arb_mul(k_x,func_eta_k->p_2,x, prec); //k*x
    arb_mul(k_y,func_eta_k->p_2,y, prec); //k*y
    arb_mul(k_eta,func_eta_k->p_2,func_eta_k->p_1, prec); //k*η
    
    arb_sqr(x_2,x,prec); //x^2
    arb_sqr(y_2,y,prec); //y^2
    
    GW_I_c_func_analyze(I_c,x,y,func_eta_k->p_2,prec); //用解析公式
    GW_I_s_func_analyze(I_s,x,y,func_eta_k->p_2,prec);
    
    
    //被积函数的具体形式见1804.07732 (43)
    //左边的xy部分
    arb_div(s,x_2,y_2,prec); // x^2/y^2 *[1- ((1+x^2-y^2)^2)/(4*x^2) ]^2
    arb_sub(t,x_2,y_2,prec);
    arb_add_ui(t,t,1,prec);
    arb_sqr(t,t,prec);
    arb_div(t,t,x_2,prec);
    arb_div_ui(t,t,4,prec);
    arb_neg(t,t);
    arb_add_ui(t,t,1,prec);
    arb_sqr(t,t,prec);
    arb_mul(s,s,t,prec);
    
    
    //中间曲率扰动功率谱部分
    arb_log(q,k_x,prec);
    power_spectrum(t, q, prec); //这里传递的是以对数的形式传入
    arb_log(q,k_y,prec);
    power_spectrum(w, q, prec);
    arb_mul(t,t,w,prec);
    arb_mul(s,s,t,prec);
    
    
    //右边括号部分
    arb_cos(t,k_eta,prec);
    arb_mul(t,t,I_c,prec);
    arb_sqr(t,t,prec);
    
    arb_sin(w,k_eta,prec);
    arb_mul(w,w,I_s,prec);
    arb_sqr(w,w,prec);
    arb_add(t,t,w,prec);
    
    arb_add(w,k_eta,k_eta,prec);
    arb_sin(w,w,prec);
    arb_mul(w,w,I_c,prec);
    arb_mul(w,w,I_s,prec);
    arb_add(t,t,w,prec);
    
    arb_mul(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(q);
    arb_clear(k_x);
    arb_clear(k_y);
    arb_clear(k_eta);
    arb_clear(x_2);
    arb_clear(y_2);
    arb_clear(I_c);
    arb_clear(I_s);
    
    return 0;
}

//诱导引力波的功率谱
int GW_power_spectra_01(arb_t res, const arb_t eta, const arb_t k, slong prec)
{
    arb_t s,t,w;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //其中η和k是积分函数的参数，需传入
    //这里，对于结构体 Func_transfer_parameter 需手动分配内存
    struct Func_transfer_parameter *func_eta_k = (struct Func_transfer_parameter *)calloc(1,sizeof(struct Func_transfer_parameter));
    
    arb_init(func_eta_k->p_1); //使用arb_t变量前初始化
    arb_init(func_eta_k->p_2);
    
    arb_set(func_eta_k->p_1,eta); //设定传递参数
    arb_set(func_eta_k->p_2,k);
    
    
    //积分前面的系数
    arb_mul(t,k,eta,prec);
    arb_sqr(t,t,prec);
    arb_mul_ui(t,t,81,prec);
    arb_inv(t,t,prec);
    arb_mul_ui(t,t,4,prec);
    
    //二元积分
    integration_binary_func(w, interior_GW_power_spectra_01_x_y, func_eta_k, 0,
                            Int_GW_power_spectra_min, Int_GW_power_spectra_max, Int_GW_power_spectra_precision, 
                            Int_GW_power_spectra_iterate_min, Int_GW_power_spectra_iterate_max,
                            Int_GW_power_spectra_01_x_y_a, NULL, 0,
                            Int_GW_power_spectra_01_x_y_b, NULL, 0,
                            Int_GW_power_spectra_precision,
                            Int_GW_power_spectra_iterate_min, Int_GW_power_spectra_iterate_max,
                            prec);
    
    arb_mul(res,t,w,prec);//乘上前面的系数
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    free(func_eta_k);
    
    return 0;
}


//当前引力波能量密度积分函数
static int interior_GW_current_energy_density_01(arb_t res, const arb_t x, const arb_t y, void* param, const slong order, slong prec)
{
    arb_t s,t,w,q,k_x,k_y,x_2,y_2,I_c,I_s;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(q);
    arb_init(k_x);
    arb_init(k_y);
    arb_init(x_2);
    arb_init(y_2);
    arb_init(I_c);
    arb_init(I_s);
    
    //接收传入参数
    struct Func_transfer_parameter *func_k; //这里不需要手动分配，只需将其指向传入的指针即可
    func_k=param;
    
    arb_mul(k_x,func_k->p_1,x, prec); //k*x
    arb_mul(k_y,func_k->p_1,y, prec); //k*y
    
    arb_sqr(x_2,x,prec); //x^2
    arb_sqr(y_2,y,prec); //y^2
    
    GW_I_c_func_analyze(I_c,x,y,func_k->p_1,prec); //用解析公式
    GW_I_s_func_analyze(I_s,x,y,func_k->p_1,prec);
    
    
    //被积函数的具体形式见1804.07732 (51)
    //左边的xy部分
    arb_div(s,x_2,y_2,prec); // x^2/y^2 *[1- ((1+x^2-y^2)^2)/(4*x^2) ]^2
    arb_sub(t,x_2,y_2,prec);
    arb_add_ui(t,t,1,prec);
    arb_sqr(t,t,prec);
    arb_div(t,t,x_2,prec);
    arb_div_ui(t,t,4,prec);
    arb_neg(t,t);
    arb_add_ui(t,t,1,prec);
    arb_sqr(t,t,prec);
    arb_mul(s,s,t,prec);
    
    
    //中间曲率扰动功率谱部分
    arb_log(q,k_x,prec);
    power_spectrum(t, q, prec); //这里传递的是以对数的形式传入
    arb_log(q,k_y,prec);
    power_spectrum(w, q, prec);
    arb_mul(t,t,w,prec);
    arb_mul(s,s,t,prec);
    
    
    //右边括号部分
    arb_sqr(w,I_c,prec);
    arb_sqr(q,I_s,prec);
    arb_add(t,w,q,prec);
    
    arb_mul(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(q);
    arb_clear(k_x);
    arb_clear(k_y);
    arb_clear(x_2);
    arb_clear(y_2);
    arb_clear(I_c);
    arb_clear(I_s);
    
    return 0;
}

//当前的引力波能量密度
int GW_current_energy_density_01(arb_t res, const arb_t k, slong prec)
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
    
    
    //积分前面的系数 c_g * Ω_{r,0} / 972
    arb_div(s,effective_g_star,effective_g_star_current,prec); //c_g= g_star/g_{star,0} * (g_{star,s,0}/g_{star,s})^{4/3}
    arb_div(t,effective_g_star_current_entropy,effective_g_star,prec);
    arb_one(w);
    arb_mul_ui(w,w,4,prec);
    arb_div_ui(w,w,3,prec);
    arb_pow(t,t,w,prec);
    arb_mul(c_g,s,t,prec);
    arb_mul(s,c_g,Omega_radiation,prec);
    arb_div_ui(s,s,972,prec);
    
    //二元积分
    integration_binary_func(w, interior_GW_current_energy_density_01, func_k, 0,
                            Int_GW_power_spectra_min, Int_GW_power_spectra_max, Int_GW_power_spectra_precision, 
                            Int_GW_power_spectra_iterate_min, Int_GW_power_spectra_iterate_max,
                            Int_GW_power_spectra_01_x_y_a, NULL, 0,
                            Int_GW_power_spectra_01_x_y_b, NULL, 0,
                            Int_GW_power_spectra_precision,
                            Int_GW_power_spectra_iterate_min, Int_GW_power_spectra_iterate_max,
                            prec);
    
    arb_mul(res,s,w,prec);//乘上前面的系数
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(c_g);
    free(func_k);
    
    return 0; 
}

