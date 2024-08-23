#include "ps_abundance_simple.h" 
#include <stdlib.h>
#include <arb_hypgeom.h>

//
//这里与前面的compaction function 方法不一样的地方在于，
//进入视界的尺度是通过质量来计算出来
//从而影响窗口函数和转移函数的值，进行影响到协方差阵Σ
//这方法可以适用于高斯和非高斯两种情况
//

static int interior_PS_abundance_simpele_compact_beta_m(arb_t res, const arb_t cl, void* p, const slong order, slong prec)
{
    arb_t t,s;
    
    arb_init(t);
    arb_init(s);
    
    
    if(Critical_Collapse_Effect==true) //考虑临界坍缩效应
    {
        //∫K*(C_l-3/8*C_l^2-C_th)^γ*P(C_l)*d(C_l)
        arb_sqr(t,cl,prec);
        arb_mul_ui(t,t,3,prec);
        arb_div_ui(t,t,8,prec);
        arb_sub(t,cl,t,prec);
        
        arb_sub(t,t,PS_C_th,prec);
        arb_abs(t,t);
        arb_pow(t,t,Mass_gamma,prec);
        
        arb_mul(t,t,Mass_K,prec);
        
        Probability_C_l(s,cl,prec); //cl对应概率 P(cl), 这里有非高斯性的影响
        
        arb_mul(res,t,s,prec);
        
    }else //不考虑临界坍缩效应
    {
        Probability_C_l(s,cl,prec); //cl对应概率 P(cl), 这里有非高斯性的影响
        
        arb_mul(res,s,Mass_K,prec);
    }
    
    
    arb_clear(s);
    arb_clear(t);
    
    return 0;
}


//利用密度扰动 δ 估算，高斯情况
int PS_abundance_simpele_compact_beta_m(arb_t res, const arb_t m, slong prec)
{
    arb_t s,R,b,sigma_xx,sigma_xy,sigma_yx,sigma_yy;
    
    arb_init(s);
    arb_init(R);
    arb_init(b);
    arb_init(sigma_xx);
    arb_init(sigma_xy);
    arb_init(sigma_yx);
    arb_init(sigma_yy);
    
    PS_variance_help_M_to_R(R,m,prec); //质量 m 对应的视界尺度 R
    
    //利用R获取协方差阵的值
    Variance_XX(sigma_xx,R,prec);
    Variance_XY(sigma_xy,R,prec);
    arb_set(sigma_yx,sigma_xy);
    Variance_YY(sigma_yy,R,prec);
    
    //通过临时改变全局变量协方差阵的方法，利用前面已有的函数
    arb_swap(sigma_xx,PS_Sigma_XX);
    arb_swap(sigma_xy,PS_Sigma_XY);
    arb_swap(sigma_yx,PS_Sigma_YX);
    arb_swap(sigma_yy,PS_Sigma_YY);
    
    arb_init(b); //积分上界 4/3
    arb_one(b);
    arb_mul_ui(b,b,4,prec);
    arb_div_ui(b,b,3,prec);
    
    Integration_arb(s, interior_PS_abundance_simpele_compact_beta_m, NULL, 0, 
                    PS_C_l_th, b, PS_abundance_simple_int_precision,
                    Integration_iterate_min, Integration_iterate_max, prec);
    arb_set(res,s);
    
    arb_swap(sigma_xx,PS_Sigma_XX);//恢复全局变量的值
    arb_swap(sigma_xy,PS_Sigma_XY);
    arb_swap(sigma_yx,PS_Sigma_YX);
    arb_swap(sigma_yy,PS_Sigma_YY);
    
    arb_clear(s);
    arb_clear(R);
    arb_clear(b);
    arb_clear(sigma_xx);
    arb_clear(sigma_xy);
    arb_clear(sigma_yx);
    arb_clear(sigma_yy);
    
    return 0;
}


static int interior_PS_abundance_simpele_compact_beta_all(arb_t res, const arb_t m, void* p, const slong order, slong prec)
{
    arb_t t,s;
    
    arb_init(t);
    arb_init(s);
    
    //以m作为积分变量
    //积分∫β(m)dln(m)=∫β(m)/m*dm
    //以m作自变量，积分函数为 β(m)/m
    
    PS_abundance_simpele_compact_beta_m(t,m,prec);
    arb_div(res,t,m,prec);
    
    arb_clear(t);
    arb_clear(s);
    
    return 0;
}

int PS_abundance_simpele_compact_beta_all(arb_t res, slong prec)
{
    arb_t t;
    arb_init(t);
    
    Integration_arb(t, interior_PS_abundance_simpele_compact_beta_all, NULL, 0, 
                    PS_abundance_simple_int_min, PS_abundance_simple_int_max, PS_abundance_simple_int_precision,
                    Integration_iterate_min, Integration_iterate_max, prec);
    arb_set(res,t);
    
    arb_clear(t);
    return 0; 
}


int PS_abundance_simpele_compact_f_m(arb_t res, const arb_t m, slong prec)
{
    arb_t beta_m,coef;
    
    arb_init(beta_m);
    arb_init(coef);
    
    beta_m_to_f_m_coefficient(coef, prec); //β 到 f 间的转换系数
    
    PS_abundance_simpele_compact_beta_m(beta_m,m,prec); //求出β(m)
    
    arb_mul(res,coef,beta_m,prec);
    
    arb_clear(beta_m);
    arb_clear(coef);
    
    return 0;
}

static int interior_PS_abundance_simpele_compact_f_all(arb_t res, const arb_t m, void* p, const slong order, slong prec)
{
    arb_t t,s;
    
    arb_init(t);
    arb_init(s);
    
    //以m作为积分变量
    //积分∫f(m)dln(m)=∫f(m)/m*dm
    //以m作自变量，积分函数为 f(m)/m
    
    PS_abundance_simpele_compact_f_m(t,m,prec);
    arb_div(res,t,m,prec);
    
    arb_clear(t);
    arb_clear(s);
    
    return 0;
}


int PS_abundance_simpele_compact_f_all(arb_t res, slong prec)
{
    arb_t t;
    arb_init(t);
    
    Integration_arb(t, interior_PS_abundance_simpele_compact_f_all, NULL, 0, 
                    PS_abundance_simple_int_min, PS_abundance_simple_int_max, PS_abundance_simple_int_precision,
                    Integration_iterate_min, Integration_iterate_max, prec);
    arb_set(res,t);
    
    arb_clear(t);
    return 0; 
}

