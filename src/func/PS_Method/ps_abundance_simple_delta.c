#include "ps_abundance_simple.h" 
#include <stdlib.h>
#include <arb_hypgeom.h>

static int interior_PS_abundance_simpele_delta_beta_m(arb_t res, const arb_t delta, void* sigma_square, const slong order, slong prec)
{
    arb_t t,s,w;
    
    arb_init(t);
    arb_init(s);
    arb_init(w);
    
    arb_mul_ui(w,sigma_square,2,prec); //2*σ^2
    
    // ∫2/sqrt(2*π*σ^2)*Exp( - δ^2/(2*σ^2) )*K*(δ-δ_c)^γ dδ
    
    arb_sqr(s,delta,prec);
    arb_div(s,s,w,prec);
    arb_neg(s,s);
    arb_exp(s,s,prec);
    
    arb_mul_ui(s,s,2,prec);
    
    arb_mul(t,w,Pi,prec);
    arb_sqrt(t,t,prec);
    
    arb_div(s,s,t,prec);
    
    arb_mul(s,s,Mass_K,prec);
    
    arb_sub(t,delta,PS_delta_th,prec);
    arb_abs(t,t);
    arb_pow(t,t,Mass_gamma,prec);
    
    arb_mul(res,s,t,prec);
    
    arb_clear(t);
    arb_clear(s);
    arb_clear(w);
    
    return 0;
}


//利用密度扰动 δ 估算，高斯情况
int PS_abundance_simpele_delta_beta_m(arb_t res, const arb_t m, slong prec)
{
    arb_t t,s,R,sigma_square;
    
    arb_init(t);
    arb_init(s);
    arb_init(R);
    arb_init(sigma_square);
    
    PS_variance_help_M_to_R(R,m,prec); //质量 m 对应的视界尺度 R
    
    //得到尺度 R 后，通过窗口函数算出对应质量的方差σ
    PS_variance_with_window_func(sigma_square, R,
                                 PS_simple_window_func, PS_simple_window_func, 1, prec); //参数 1 功率谱使用 P_δ(k)
    
    if(Critical_Collapse_Effect==true)
    {
        //考虑临界效应，需做一个积分
        arb_t b; //积上界
        arb_init(b);
        
        arb_mul_ui(b,PS_delta_th,5,prec); //5*δ_c
        
        //需传入参数σ^2
        Integration_arb(s, interior_PS_abundance_simpele_delta_beta_m, sigma_square, 0, 
                        PS_delta_th, b, PS_abundance_simple_int_precision,
                        Integration_iterate_min, Integration_iterate_max, prec);
        
        arb_clear(b);
    }else
    {
        //不考虑临界效应时
        //再通过高斯分布积分（互补误差函数）求出对应的概率
        //β(M)=K*erfc(δ_c/sqrt(2*σ^2))
        
        arb_mul_ui(t,sigma_square,2,prec);
        arb_sqrt(t,t,prec);
        arb_div(t,PS_delta_th,t,prec);
        
        arb_hypgeom_erfc(s,t,prec/2); //降点精度，加快速度
        
        arb_mul(s,s,Mass_K,prec); 
        
    }
    
    arb_set(res,s);
    
    arb_clear(t);
    arb_clear(s);
    arb_clear(R);
    arb_clear(sigma_square);
    
    return 0;
}


static int interior_PS_abundance_simpele_delta_beta_all(arb_t res, const arb_t m, void* p, const slong order, slong prec)
{
    arb_t t,s;
    
    arb_init(t);
    arb_init(s);
    
    //以m作为积分变量
    //积分∫β(m)dln(m)=∫β(m)/m*dm
    //以m作自变量，积分函数为 β(m)/m
    
    PS_abundance_simpele_delta_beta_m(t,m,prec);
    arb_div(res,t,m,prec);
    
    arb_clear(t);
    arb_clear(s);
    
    return 0;
}

int PS_abundance_simpele_delta_beta_all(arb_t res, slong prec)
{
    arb_t t;
    arb_init(t);
    
    Integration_arb(t, interior_PS_abundance_simpele_delta_beta_all, NULL, 0, 
                    PS_abundance_simple_int_min, PS_abundance_simple_int_max, PS_abundance_simple_int_precision,
                    Integration_iterate_min, Integration_iterate_max, prec);
    arb_set(res,t);
    
    arb_clear(t);
    return 0; 
}


int PS_abundance_simpele_delta_f_m(arb_t res, const arb_t m, slong prec)
{
    arb_t beta_m,coef;
    
    arb_init(beta_m);
    arb_init(coef);
    
    beta_m_to_f_m_coefficient(coef, prec); //β 到 f 间的转换系数
    
    PS_abundance_simpele_delta_beta_m(beta_m,m,prec); //求出β(m)
    
    arb_mul(res,coef,beta_m,prec);
    
    arb_clear(beta_m);
    arb_clear(coef);
    
    return 0;
}

static int interior_PS_abundance_simpele_delta_f_all(arb_t res, const arb_t m, void* p, const slong order, slong prec)
{
    arb_t t,s;
    
    arb_init(t);
    arb_init(s);
    
    //以m作为积分变量
    //积分∫f(m)dln(m)=∫f(m)/m*dm
    //以m作自变量，积分函数为 f(m)/m
    
    PS_abundance_simpele_delta_f_m(t,m,prec);
    arb_div(res,t,m,prec);
    
    arb_clear(t);
    arb_clear(s);
    
    return 0;
}


int PS_abundance_simpele_delta_f_all(arb_t res, slong prec)
{
    arb_t t;
    arb_init(t);
    
    Integration_arb(t, interior_PS_abundance_simpele_delta_f_all, NULL, 0, 
                    PS_abundance_simple_int_min, PS_abundance_simple_int_max, PS_abundance_simple_int_precision,
                    Integration_iterate_min, Integration_iterate_max, prec);
    arb_set(res,t);
    
    arb_clear(t);
    return 0; 
}

