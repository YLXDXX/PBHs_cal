#include "ps_abundance.h" 
#include <stdlib.h>


//在PBHs生成时，PBHs的能量密度占总能量密度的比记为 β
//下面求 β 与 m 间的关系 β(m)


static int interior_PS_abundance_beta_m(arb_t res, const arb_t cl, void* p, const slong order, slong prec)
{
    //不考虑临界坍缩效应时用，这里积分后，求总β便不用再积分
    arb_t t;
    arb_init(t);
    
    Probability_C_l(t,cl,prec); //cl对应概率 P(cl)
    
    arb_set(res,t);
    
    arb_clear(t);
    
    return 0;
}


int PS_abundance_beta_m(arb_t res, const arb_t m, slong prec)
{
    arb_t t,s,w,cl;
    
    arb_init(t);
    arb_init(s);
    arb_init(w);
    arb_init(cl);
    
    if(Critical_Collapse_Effect==true) //考虑临界坍缩效应
    {
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
        
    }else //不考虑临界坍缩效应
    {
        //此时质量 m 仅能取一个值，M=K·M_H, 考虑相对质量，m的取值为 m=K
        if( arb_eq(m,Mass_K) )
        {
            //β=β(M)=K·∫P(C_l)dC_l
            
            arb_t a,b;
            arb_init(a);
            arb_init(b);
            
            arb_set(a,PS_C_l_th); //积分下界 C_{l,th}
            arb_one(b); //积分上界 4/3
            arb_mul_ui(b,b,4,prec);
            arb_div_ui(b,b,3,prec);
            
            Integration_arb(t, interior_PS_abundance_beta_m, NULL, 0, 
                            a, b, PS_abundance_int_precision,
                            Integration_iterate_min,Integration_iterate_max, prec);
            arb_set(res,t);
            
            arb_clear(a);
            arb_clear(b);
        }else
        {
            printf("不考虑临界坍缩效应，M=K·M_H，请设相对质量为 K\n ");
            exit(1);
        }
    }
    
    
    arb_clear(t);
    arb_clear(s);
    arb_clear(w);
    arb_clear(cl);
    
    return 0;
}

static int interior_PS_abundance_beta_all(arb_t res, const arb_t m, void* p, const slong order, slong prec)
{
    arb_t t,s;
    
    arb_init(t);
    arb_init(s);
    
    //以m作为积分变量
    //积分∫β(m)dln(m)=∫β(m)/m*dm
    //以m作自变量，积分函数为 β(m)/m
    
    PS_abundance_beta_m(t,m,prec);
    
    arb_div(res,t,m,prec);
    
    
    arb_clear(t);
    arb_clear(s);
    
    return 0;
}



int PS_abundance_beta_all(arb_t res, slong prec)
{
    arb_t t,a,b;
    
    arb_init(t);
    arb_init(a);
    arb_init(b);
    
    if(Critical_Collapse_Effect==true) //考虑临界坍缩效应
    {
        //∫_0^{M_ratio_max} f(m)dln(m)
        //积分下界不能取到0，有除以零的运算
        arb_set(a,PS_M_ratio_min);
        arb_set(b,PS_M_ratio_max);
        
        //使用新的gauss_kronrod积分算法
        Integration_arb(t, interior_PS_abundance_beta_all, NULL, 0, 
                        a, b, PS_abundance_int_precision,
                        Integration_iterate_min,Integration_iterate_max, prec);
        arb_set(res,t);
        
    }else //不考虑临界坍缩效应
    {
        //β=β(M)=K·∫P(C_l)dC_l
        arb_set(t,Mass_K);
        PS_abundance_beta_m(res,t,prec);
    }
    
    
    
    arb_clear(t);
    arb_clear(a);
    arb_clear(b);
    
    return 0;
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



static int interior_PS_abundance_f_all(arb_t res, const arb_t m, void* p, const slong order, slong prec)
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
    arb_t t,a,b;
    
    arb_init(t);
    arb_init(a);
    arb_init(b);
    
    if(Critical_Collapse_Effect==true) //考虑临界坍缩效应
    {
        //∫_0^{M_ratio_max} f(m)dln(m)
        //积分下界不能取到0，有除以零的运算
        arb_set(a,PS_M_ratio_min);
        arb_set(b,PS_M_ratio_max);
        
        //使用新的gauss_kronrod积分算法
        Integration_arb(t, interior_PS_abundance_f_all, NULL, 0, 
                        a, b, PS_abundance_int_precision,
                        Integration_iterate_min,Integration_iterate_max, prec);
        arb_set(res,t);
        
    }else //不考虑临界坍缩效应
    {
        //f=f(M)
        arb_set(t,Mass_K);
        PS_abundance_f_m(res,t,prec);
    }
    
    
    arb_clear(t);
    arb_clear(a);
    arb_clear(b);
    
    return 0;
}

