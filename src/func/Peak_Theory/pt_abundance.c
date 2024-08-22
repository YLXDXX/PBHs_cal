#include "pt_abundance.h"
#include <stdlib.h>

//当有了数密度后，β 或 f 都可以直接得到
//这里需要区分的是：是否使用相对质量
//注意，若直接使用数密度得到 f 则不能使用相对质量
//另外，这里的相对质量跟PS方法的相对质量不太一样
//PS：M/M_H, PT: M/M(k_*)

static void into_horizon_V(arb_t res, const arb_t r_m, slong prec)
{
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    //采用最简单的，4/3*π*r_m^3
    
    arb_one(s);
    arb_mul_ui(s,s,4,prec);
    arb_div_ui(s,s,3,prec);
    arb_mul(s,s,Pi,prec);
    arb_pow_ui(t,r_m,3,prec);
    
    arb_mul(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
}

//PBH abundance β_PBH(M)
int PT_abundance_beta_m(arb_t res, const arb_t M, slong prec)
{
    arb_t s,t,w;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    if(PT_Mass_Relative==true) //使用相对质量
    {
        //β(M/M_H) dln(M/M_H) = M/M_H * n * V(r_m) dln(M/M_H)
        // M/M_* --> M
        PBH_number_density_M(s,M,prec); //此函数可以直接处理相对质量和非相对质量的M
        
        into_horizon_V(t,R_MAX,prec);
        arb_mul(s,s,t,prec);
        
        //临界坍缩效应 M=M_H*K*(C-C_th)^γ=M_* * (r_m*k_star)^2 * exp(2*ζ_m) * K*(C-C_th)^γ
        // M/M_H=M/M_star * 1/ [ (r_m*k_star)^3 * exp(2*ζ_m)]
        // β(M/M_*) dln(M/M_*) = M/M_star *  1/ [ (r_m*k_star)^2 * exp(2*ζ_m)] * n * V(r_m) dln(M/M_H)
        
        // M/M_star/[ (r_m*k_star)^2 * exp(2*ζ_m)]
        arb_mul(t,R_MAX,K_star,prec);//相应的r_m在调用PBH_number_density_M函数时已求出
        arb_sqr(t,t,prec);
        arb_mul(t,t,PT_beta_cal_need_exp_2_zeta_m,prec);
        arb_div(t,M,t,prec);
        
        arb_mul(res,s,t,prec); //乘以相对质量
        
        
    }else //不使用相对质量
    {
        //β(M) dln(M) = M/M_H * n * V(r_m) dln(M)
        
        PBH_number_density_M(s,M,prec); //此函数可以直接处理相对质量和非相对质量的M
        
        //M_H= (r_m*k_star)^2 * exp(2*ζ_m) * M_star
        arb_mul(t,R_MAX,K_star,prec);//相应的r_m在调用PBH_number_density_M函数时已求出
        arb_sqr(t,t,prec);
        arb_mul(t,t,PT_beta_cal_need_exp_2_zeta_m,prec);
        Horizon_reentry_k_to_M(w,K_star,prec);
        arb_mul(t,t,w,prec);
        
        arb_div(t,M,t,prec);
        into_horizon_V(w,R_MAX,prec);
        
        arb_mul(s,s,t,prec);
        arb_mul(res,s,w,prec);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    
    return 0;
}

static int interior_PT_abundance_beta_all(arb_t res, const arb_t m, void* p, const slong order, slong prec)
{
    arb_t t;
    arb_init(t);
    
    //注意，这里是否使用还对质量，积分表达式都是一样的
    //否使用还对质量的区别在于积分的上下限
    
    //以m作为积分变量
    //积分∫β(m)dln(m)=∫β(m)/m*dm
    //以m作自变量，积分函数为 β(m)/m
    
    PT_abundance_beta_m(t,m,prec);
    
    arb_div(res,t,m,prec);
    
    arb_clear(t);
    
    return 0;
}

int PT_abundance_beta_all(arb_t res, slong prec)
{
    arb_t s,t,a,b;
    arb_init(s);
    arb_init(t);
    arb_init(a);
    arb_init(b);
    
    //否使用还对质量的区别在于积分的上下限
    if(PT_Mass_Relative==true) //使用相对质量
    {
        //∫_0^{M_ratio_max} f(m)dln(m)
        //积分下界不能取到0，有除以零的运算
        arb_set(a,PT_M_ratio_min);
        arb_set(b,PT_M_ratio_max); //由于peak theory与PS方法相对质量利用的质量不同，其取值范围也不同
        //arb_set_str(b,"1.6",prec);
        
        int int_n_min,int_n_max;
        if(Integral_method==double_exponential)
        {
            int_n_min=4;
            int_n_max=6;
            
        }else if(Integral_method==gauss_kronrod_iterate)
        {
            int_n_min=2;
            int_n_max=4;
        }else
        {
            printf(" Peak_Theory -> abundance.c -> PT_abundance_beta_all -> 积分迭代次数设置错误\n");
            exit(1);
        }
        
        //使用新的gauss_kronrod积分算法
        Integration_arb(res, interior_PT_abundance_beta_all, NULL, 0, 
                        a, b, PS_abundance_int_precision,
                        int_n_min,int_n_max, prec);
        
    }else //不使用相对质量
    {
        printf("不使用相对质量时，积分区间太大，请使用相对质量求 β_all\n");
        exit(1);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(a);
    arb_clear(b);
    
    return 0;
}






//PBH abundance f_PBH(M)


int PT_abundance_f_m(arb_t res, const arb_t m, slong prec)
{
    arb_t beta_m,coef;
    
    arb_init(beta_m);
    arb_init(coef);
    
    beta_m_to_f_m_coefficient(coef, prec); //β 到 f 间的转换系数
    
    //求出β(m)
    PT_abundance_beta_m(beta_m,m,prec);
    
    arb_mul(res,coef,beta_m,prec);
    
    arb_clear(beta_m);
    arb_clear(coef);
    
    return 0;
}


static int interior_PT_abundance_f_all(arb_t res, const arb_t m, void* p, const slong order, slong prec)
{
    arb_t t;
    arb_init(t);
    
    //以m作为积分变量
    //积分∫f(m)dln(m)=∫f(m)/m*dm
    //以m作自变量，积分函数为 f(m)/m
    
    PT_abundance_f_m(t,m,prec);
    
    arb_div(res,t,m,prec);
    
    
    arb_clear(t);
    
    return 0;
}

int PT_abundance_f_all(arb_t res, slong prec)
{
    arb_t s,t,a,b;
    arb_init(s);
    arb_init(t);
    arb_init(a);
    arb_init(b);
    
    //否使用还对质量的区别在于积分的上下限
    if(PT_Mass_Relative==true) //使用相对质量
    {
        //∫_0^{M_ratio_max} f(m)dln(m)
        //积分下界不能取到0，有除以零的运算
        arb_set(a,PT_M_ratio_min);
        arb_set(b,PT_M_ratio_max); //由于peak theory与PS方法相对质量利用的质量不同，其取值范围也不同
        //arb_set_str(b,"1.6",prec);
        
        int int_n_min,int_n_max;
        if(Integral_method==double_exponential)
        {
            int_n_min=4;
            int_n_max=6;
            
        }else if(Integral_method==gauss_kronrod_iterate)
        {
            int_n_min=2;
            int_n_max=4;
        }else
        {
            printf(" Peak_Theory -> abundance.c -> PT_abundance_beta_all -> 积分迭代次数设置错误\n");
            exit(1);
        }
        
        //使用新的gauss_kronrod积分算法
        Integration_arb(res, interior_PT_abundance_f_all, NULL, 0, 
                        a, b, PS_abundance_int_precision,
                        int_n_min,int_n_max, prec);
        
    }else //不使用相对质量
    {
        printf("不使用相对质量时，积分区间太大，请使用相对质量求 f_all\n");
        exit(1);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(a);
    arb_clear(b);
    
    return 0;
}





