#include "abundance.h"
#include <stdlib.h>


//PBH abundance f_PBH(M)
int PBH_abundance_f_to_M(arb_t res,const arb_t M, slong prec)
{
    //函数中所用变量
    arb_t s,t,w,mu;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(mu);
    
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
            if( Relative_Mass )
            {
                //使用相对质量
                arb_set_str(t,"1E20",prec);
                arb_div(s,M,t,prec);
                Horizon_reentry_mu_M(mu,s,prec);//通过M获得µ
            }else
            {
                //不使用相对质量
                Horizon_reentry_mu_M(mu,M,prec);//通过M获得µ
            }
            
            
            
            Horizon_reentry_D_ln_M_to_mu(s,mu,prec);
            arb_inv(s,s,prec);
            arb_mul(w,w,s,prec);
            
            arb_sqrt(s,Power_A,prec);
            arb_div(s,mu,s,prec);
            N_pk_f_xi(t,s,prec);
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


