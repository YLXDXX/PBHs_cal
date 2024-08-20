#include "ps_generate_mass.h"
#include <stdlib.h>

//黑洞质量m与C_l的关系，已知 C_l 求 m
int PS_C_l_to_M_ratio(arb_t res, const arb_t cl, slong prec)
{
    arb_t t,s;
    arb_init(t);
    arb_init(s);
    
    //这里是以质量的比率 M/M_* 作为因变量的 m=M/M_*
    //M=(r_m*k_star)^2 * exp(2*ζ_m)*K*(C_l-3/8*(C_l)^2-C_th)^γ*M_(k_star)
    //M=K*(C_l-3/8*(C_l)^2-C_th)^γ*M_H
    
    //在PS方法中，对于视界质量的计算，忽略扰动ζ的影响
    //有：m=K*(C_l-3/8*(C_l)^2-C_th)^γ
    
    arb_sqr(s,cl,prec);
    arb_mul_ui(s,s,3,prec);
    arb_div_ui(s,s,8,prec);
    arb_sub(s,cl,s,prec);
    arb_sub(s,s,PS_C_th,prec);
    
    arb_pow(t,s,Mass_gamma,prec);
    arb_mul(res,t,Mass_K,prec);
    
    /*
    zeta_profile_n(s,R_MAX,0,prec); // ζ_m
    arb_mul_ui(s,s,2,prec);
    arb_exp(s,s,prec);
    arb_mul(t,t,s,prec);
    
    arb_mul(s,K_star,R_MAX,prec);
    arb_sqr(s,s,prec);
    arb_mul(res,t,s,prec);
    */
    
    
    arb_clear(t);
    arb_clear(s);
    
    return 0;
}


//黑洞质量m与C_l的关系，已知 m 求 C_l
int PS_M_ratio_to_C_l(arb_t res, const arb_t m, slong prec)
{
    //这里是以质量的比率 M/M_H 作为自变量的 m=M/M_H
    //利用m与C_l的关系，反解C_l，为一元二次方程，取最小的根
    
    if ( arb_gt(m,PS_M_ratio_max) || arb_lt(m,PS_M_ratio_min) )
    {
        printf("PS_Method -> ps_generate_mass -> PS_M_ratio_to_C_l 超出 I 型扰动质量范围：");
        arb_printn(PS_M_ratio_min, 15,0);printf("\n");
        arb_printn(PS_M_ratio_max, 15,0);printf("\n");
        arb_printn(m, 15,0);printf("\n");
        exit(1);
    }
    
    arb_t t,s,b,c;
    
    arb_init(t);
    arb_init(s);
    arb_init(b);
    arb_init(c);
    
    
    //概据两者的函数关系，m的最大值在4/3处取得
    //有两种表述
    // C^2 - 8/3*C + 8/3*[ [m*(k_star*r_m)^{-2}*e^(-2*ζ_m)*K^{-1}]^{1/Gamma} +C_th ]
    // C^2 - 8/3*C + 8/3*[ (m/K)^{1/Gamma} +C_th ] //现使用下面使用此种表述
    
    //二次函数的系数 b,c
    arb_one(b);
    arb_mul_ui(b,b,8,prec);
    arb_div_ui(b,b,3,prec);
    arb_neg(b,b);
    
    /*
    arb_mul(s,K_star,R_MAX,prec); //此处的 R_MAX 需要预先求出
    arb_sqr(s,s,prec);
    arb_inv(s,s,prec);
    arb_mul(s,s,m,prec); //m*(k_star*r_m)^{-2}
    
    zeta_profile_n(t,R_MAX,0,prec); // ζ_m
    arb_mul_ui(t,t,2,prec);
    arb_neg(t,t);
    arb_exp(t,t,prec);
    arb_div(t,t,Mass_K,prec);
    arb_mul(t,t,s,prec); //m*(k_star*r_m)^{-2}*e^(-2*ζ_m)*K^{-1}
    
    arb_inv(s,Mass_gamma,prec);
    arb_pow(c,t,s,prec);
    arb_add(c,c,PS_C_th,prec);
    arb_mul_ui(c,c,8,prec);
    arb_div_ui(c,c,3,prec);
    */
    
    // c采用 8/3*[ (m/K)^{1/Gamma} +C_th ]
    
    arb_div(t,m,Mass_K,prec);
    arb_inv(s,Mass_gamma,prec);
    arb_pow(t,t,s,prec);
    arb_add(c,t,PS_C_th,prec);
    arb_mul_ui(c,c,8,prec);
    arb_div_ui(c,c,3,prec);
    
    
    //求出b,c后，用求根公式算
    arb_mul_si(s,c,4,prec); //a=1
    arb_neg(s,s);
    arb_sqr(t,b,prec);
    arb_add(t,t,s,prec); 
    arb_sqrt(t,t,prec); // sqrt(b^2-4ac)
    
    arb_neg(s,b);
    arb_sub(s,s,t,prec);//-b - sqrt(b^2-4ac)
    
    arb_div_si(res,s,2,prec);
    
    
    arb_clear(t);
    arb_clear(s);
    arb_clear(b);
    arb_clear(c);
    
    return 0;
}


//获得相对质量取值的范围
void PS_get_relative_M_range(arb_t min, arb_t max, slong prec)
{
    arb_t t,s;
    arb_init(t);
    arb_init(s);
    
    
    if(Critical_Collapse_Effect==true) //考虑临界坍缩效应
    {
        //计算相对质量的最大值
        arb_one(t); // t = 4/3
        arb_mul_ui(t,t,4,prec);
        arb_div_ui(t,t,3,prec); // t = 4/3 ，对应于两类型PBH的分界点
        PS_C_l_to_M_ratio(max,t,prec);
        
        
        //注意，这里M/M_{H}不能精确的取到最大值，会有 P(C_l->4/3) --> P(C) 的发散问题
        //但P(C)的积分没有发散问题，为了数值计算上的方便，做一点微小的截断
        arb_set_str(t, "1E-45", prec); // C_l_max=C_l_max-10^{-40}
        arb_sub(max,max,t,prec);
        
        arb_zero(min); //相对质量最小值为零
        arb_add(min,min,t,prec); //一点微小截断，避免可能的发散
        
    }else //不考虑临界坍缩效应
    {
        //M=K·M_H ⇒ min=max=K
        arb_set(min,Mass_K);
        arb_set(max,min);
    }
    
    arb_clear(t);
    arb_clear(s);
}

