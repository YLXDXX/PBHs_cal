#include "abundance.h"
#include <stdlib.h>

//当有了数密度后，β 或 f 都可以直接得到
//这里需要区分的是：是否使用相对质量
//注意，若直接使用数密度得到 f 则不能使用相对质量


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
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    if(PT_Mass_Relative==true) //使用相对质量
    {
        //β(M/M_H) dln(M/M_H) = M/M_H * n * V(r_m) dln(M/M_H)
        // M/M_H --> M
        PBH_number_density_M(s,M,prec);
        
        //注意到，这里生成质量的临界坍缩效应，在计算数密度的时候，就已经考虑过了
        into_horizon_V(t,R_MAX,prec);
        arb_mul(s,s,t,prec);
        
        arb_mul(res,s,M,prec);
        
    }else //不使用相对质量
    {
        
    }
    
    arb_clear(s);
    arb_clear(t);
    
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
        arb_set_str(a,"1E-5",prec); //a不能精确等于0
        arb_set(b,PS_M_ratio_max);
        //arb_set_str(b,"1.6",prec);
        
        //使用新的gauss_kronrod积分算法
        Integration_arb(res, interior_PT_abundance_beta_all, NULL, 0, 
                        a, b, PS_abundance_f_all_precision,
                        Integration_iterate_min,Integration_iterate_max, prec);
        
    }else //不使用相对质量
    {
        
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(a);
    arb_clear(b);
    
    return 0;
}






//PBH abundance f_PBH(M)

int PT_abundance_f_all(arb_t res, slong prec)
{
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    if(PT_Mass_Relative==true) //使用相对质量
    {
        
    }else //不使用相对质量
    {
        
    }
    
    arb_clear(s);
    arb_clear(t);
    
    return 0;
}


int PT_abundance_f_m(arb_t res, const arb_t M, slong prec)
{
    //函数中所用变量
    arb_t s,t,w,mu,ttttttt;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(mu);
    arb_init(ttttttt);
    
    //跟功率谱有关，对于delta型式的功率谱，可解析求解
    //功率谱类型判断
    switch(Power_spectrum_type) 
    {
        case lognormal_type :
            
            //先简要计算一下
            // f_PBH(M) = M*n_PBH(M) / A 其中 A 为常数 
            //先不考虑  A 
            PBH_number_density_M(s,M,prec);
            
            arb_mul(res,M,s,prec);
            
            break;
        case delta_type :
            //在delta谱的情况下，f_PBH(M)可解析求解
            //前面常系数
            arb_sqr(s,Hubble_constant_h,prec);
            arb_mul(s,s,Omega_DM,prec);
            arb_set_str(t,"0.12",prec);
            arb_div(w,t,s,prec);
            
            //关于M
            arb_set_str(t,"1E20",prec);
            arb_div(s,M,t,prec);
            arb_mul(w,w,s,prec);
            
            //关于k_*
            arb_set_str(t,"1.56E13",prec);
            arb_div(s,K_star,t,prec);
            arb_pow_ui(t,s,3,prec); //三次方
            arb_mul(w,w,t,prec);
            
            //最后
            
            //这里需要判断一下，是否使用相对质量
            if( PT_Mass_Relative )
            {
                //使用相对质量
                arb_set_str(t,"1E20",prec);
                arb_div(s,M,t,prec);
                Horizon_reentry_M_to_mu(mu,s,ttttttt,prec);//通过M获得µ
            }else
            {
                //不使用相对质量
                Horizon_reentry_M_to_mu(mu,M,ttttttt,prec);//通过M获得µ
            }
            
            
            
            Horizon_reentry_derivative_ln_M_mu(s,mu,ttttttt,prec);
            arb_inv(s,s,prec);
            arb_mul(w,w,s,prec);
            
            arb_sqrt(s,Power_A,prec);
            arb_div(s,mu,s,prec);
            N_pk_help_f_xi(t,s,prec);
            arb_mul(w,w,t,prec);
            
            //最后P_G
            //系数
            arb_mul(s,Pi_2,Power_A,prec);
            arb_sqrt(s,s,prec);
            arb_inv(s,s,prec);
            arb_mul(w,w,s,prec);
            
            //指数
            arb_mul_si(t,Power_A,2,prec);
            arb_sqr(s,mu,prec);
            arb_div(t,s,t,prec);
            arb_neg(t,t);
            arb_exp(t,t,prec);
            arb_mul(w,w,t,prec);
            
            arb_set_str(t,"5.3E-16",prec);
            arb_div(res,w,t,prec);
            
            break;
            
            default :
                printf("Peak_Theory -> abundance -> power_spectrum_type 有误\n");
                exit(1);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(mu);
    
    return 0;
}


