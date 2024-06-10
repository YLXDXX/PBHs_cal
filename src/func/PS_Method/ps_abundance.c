#include "ps_abundance.h" 
#include <stdlib.h>


//在PBHs生成时，PBHs的能量密度占总能量密度的比记为 β
//下面求 β 与 m 间的关系 β(m)
int PS_abundance_beta_m(arb_t res, const arb_t m, slong prec)
{
    arb_t t,s,w,cl;
    
    arb_init(t);
    arb_init(s);
    arb_init(w);
    arb_init(cl);
    
    //β(m)=m*|dln(m)/dC_l|^{-1} * P(C_l) 
    //    =K*[C_l-3/8*C_l^2-C_th]^{γ+1} /[γ*(1-3/4*C_l)] * P(C_l)
    
    //首先由 m=M/M_* 解出 C_l
    PS_M_ratio_to_C_l(cl,m,prec); //由m反解出cl
    
    //分母 γ*(1-3/4*C_l)
    arb_mul_ui(s,cl,3,prec);
    arb_div_ui(s,s,4,prec);
    arb_neg(s,s);
    arb_add_si(s,s,1,prec);
    arb_mul(s,s,Mass_gamma,prec);
    
    //分子[C_l-3/8*C_l^2-C_th]^{γ+1}
    arb_sqr(t,cl,prec);
    arb_mul_ui(t,t,3,prec);
    arb_div_ui(t,t,8,prec);
    arb_sub(t,cl,t,prec);
    arb_sub(t,t,PS_C_th,prec);
    if(arb_is_negative(t))
    {
        printf("PS_Method -> ps_abundance -> C_l 小于阈值 C_l,th，请重新检查计算\n");
        exit(1);
    }
    
    arb_add_ui(w,Mass_gamma,1,prec);
    arb_pow(t,t,w,prec);
    
    arb_div(t,t,s,prec);//分子除以分母
    
    arb_mul(t,t,Mass_K,prec);//最前面系数K
    
    Probability_C_l(s,cl,prec); //cl对应概率 P(cl)
    
    arb_mul(res,t,s,prec);
    
    
    arb_clear(t);
    arb_clear(s);
    arb_clear(w);
    arb_clear(cl);
    
    return 0;
}

int interior_PS_abundance_beta_all(arb_t res, const arb_t m, void* p, const slong order, slong prec)
{
    arb_t t,s;
    
    arb_init(t);
    arb_init(s);
    
    //以m作为积分变量
    //积分∫β(m)dln(m)=∫β(m)/m*dm
    //以m作自变量，积分函数为 f(m)/m
    
    PS_abundance_beta_m(t,m,prec);
    
    arb_div(res,t,m,prec);
    
    
    arb_clear(t);
    arb_clear(s);
    
    return 0;
}



int PS_abundance_beta_all(arb_t res, slong prec)
{
    arb_t t,s,a,b;
    
    arb_init(t);
    arb_init(s);
    arb_init(a);
    arb_init(b);
    
    
    //∫_0^{M_ratio_max} f(m)dln(m)
    //积分下界不能取到0，有除以零的运算
    arb_set_str(a,"1E-90",prec); //a不能精确等于0
    arb_set(b,PS_M_ratio_max);
    
    //使用新的gauss_kronrod积分算法
    Integration_arb(res, interior_PS_abundance_beta_all, NULL, 0, 
                              a, b, PS_abundance_f_all_precision,
                              Integration_iterate_min,Integration_iterate_max, prec);
    
    
    arb_clear(t);
    arb_clear(s);
    arb_clear(a);
    arb_clear(b);
    
    return 0;
}


//生成时的能量密度 β(m)/β 和当今的能量密度 f(m)/f 间的转换系数
void beta_m_to_f_m_coefficient(arb_t res, slong prec)
{
    arb_t t,s,w;
    arb_init(t);
    arb_init(s);
    arb_init(w);
    
    switch (beta_to_f_type)
    {
        case beta_f_general_I :
            //参与 2109.00791 (A.15) 第一行
            //f(m) = Ω_m/Ω_DM * T_{r_m}/T_eq * β(M)
            
            //r_m对应尺度 k_m,需要求出 k_m 对应的温度即可
            //关于温度的k的转换关系 见 1812.00674 (7)
            // k_m/k_eq = 2(sqrt{2}-1)*...*T/T_eq
            // k_m/k_eq=1/(r_m*k_eq)
            // T/T_eq = (g_{*,eq}/g_*)^{1/2} * (g_{*,s}/g_{*,s,eq})^{1/3} * 1/( 2*(sqrt{2}-1)*r_m*k_eq )
            arb_one(w); //(g_{*,eq}/g_*)^{1/2}
            arb_div_ui(w,w,2,prec);
            arb_div(s,effective_g_star_eq,effective_g_star,prec);
            arb_pow(s,s,w,prec);
            
            
            arb_one(w); //(g_{*,s}/g_{*,s,eq})^{1/3}
            arb_div_ui(w,w,3,prec);
            arb_div(t,effective_g_star_entropy,effective_g_star_eq_entropy,prec);
            arb_pow(t,t,w,prec);
            arb_mul(s,s,t,prec);
            
            arb_one(w); //2*(sqrt{2}-1)*r_m*k_eq
            arb_mul_ui(w,w,2,prec);
            arb_sqrt(w,w,prec);
            arb_sub_ui(w,w,1,prec);
            arb_mul_ui(w,w,2,prec);
            arb_mul(w,w,R_MAX,prec);
            arb_mul(w,w,K_scale_eq,prec);
            
            arb_div(res,s,w,prec);
            
            break;
        case beta_f_general_II :
            //参与 2109.00791 (A.15) 第二行
            arb_sqr(t,Hubble_constant_h,prec);
            arb_mul(t,t,Omega_DM,prec);
            arb_set_str(w,"0.12",prec);
            arb_div(t,w,t,prec);
            
            arb_set_str(w,"106.75",prec);
            arb_div(s,effective_g_star,w,prec);
            arb_one(w);
            arb_div_ui(w,w,6,prec);
            arb_pow(s,s,w,prec);
            arb_mul(t,t,s,prec);
            
            arb_set_str(w,"2.74",prec);
            arb_div(s,w,R_m_times_K,prec);
            arb_mul(t,t,s,prec);
            
            arb_set_str(w,"1.56E13",prec);
            arb_div(s,K_star,w,prec);
            arb_mul(t,t,s,prec);
            
            arb_set_str(w,"6.88E-16",prec);
            
            arb_div(res,t,w,prec);
            
            break;
        case beta_f_myself :
            //自己推导的相关转换系数公式
            
            //f(M)=Ω_m/Ω_DM * (g_*/g_{*,eq})^(1/3) * (k_*/k_eq) * 1/[2*(√2-1)]^4 * β(M)
            arb_div(t,Omega_M,Omega_DM,prec);
            arb_div(s,effective_g_star,effective_g_star_eq,prec);
            arb_one(w);
            arb_div_ui(w,w,3,prec);
            arb_pow(s,s,w,prec);
            arb_mul(t,t,s,prec);
            
            arb_div(s,K_star,K_scale_eq,prec); //k_*/k_eq
            arb_mul(t,t,s,prec);
            
            arb_one(s); // 1/[2*(√2-1)]^4
            arb_mul_ui(s,s,2,prec);
            arb_sqrt(s,s,prec);
            arb_sub_ui(s,s,1,prec);
            arb_mul_ui(s,s,2,prec);
            arb_pow_ui(s,s,4,prec);
            
            arb_div(res,t,s,prec);
            
            break;
        default :
            printf("beta_to_f_type 输入有误\n");
            exit(1);
    }
    
    arb_clear(t);
    arb_clear(s);
    arb_clear(w);
}


//在当今时刻，PBHs的能量密度占暗能量密度的比记为 f
//下面求 f 与 m 间的关系 f(m)
int PS_abundance_f_m(arb_t res, const arb_t m, slong prec)
{
    arb_t beta_m,coef;
    
    arb_init(beta_m);
    arb_init(coef);
    
    beta_m_to_f_m_coefficient(coef, prec); //β 到 f 间的转换系数
    
    //求出β(m)
    PS_abundance_beta_m(beta_m,m,prec);
    
    arb_mul(res,coef,beta_m,prec);
    
    arb_clear(beta_m);
    arb_clear(coef);
    return 0;
}



int interior_PS_abundance_f_all(arb_t res, const arb_t m, void* p, const slong order, slong prec)
{
    arb_t t,s;
    
    arb_init(t);
    arb_init(s);
    
    //以m作为积分变量
    //积分∫f(m)dln(m)=∫f(m)/m*dm
    //以m作自变量，积分函数为 f(m)/m
    
    PS_abundance_f_m(t,m,prec);
    
    arb_div(res,t,m,prec);
    
    
    arb_clear(t);
    arb_clear(s);
    
    return 0;
}


//在当今时刻，PBHs的能量密度占暗能量密度的比记为 f
//下面求 f all
int PS_abundance_f_all(arb_t res, slong prec)
{
    arb_t t,s,a,b;
    
    arb_init(t);
    arb_init(s);
    arb_init(a);
    arb_init(b);
    
    
    //∫_0^{M_ratio_max} f(m)dln(m)
    //积分下界不能取到0，有除以零的运算
    arb_set_str(a,"1E-90",prec); //a不能精确等于0
    arb_set(b,PS_M_ratio_max);
    
    //使用新的gauss_kronrod积分算法
    Integration_arb(res, interior_PS_abundance_f_all, NULL, 0, 
                              a, b, PS_abundance_f_all_precision,
                              Integration_iterate_min,Integration_iterate_max, prec);
    
    
    arb_clear(t);
    arb_clear(s);
    arb_clear(a);
    arb_clear(b);
    
    return 0;
}

