#include "number_density.h"
#include <stdlib.h>
#include <arb_hypgeom.h> //误差函数用

//峰数密度相关函数 f(ξ)
int N_pk_f_xi(arb_t res, const arb_t xi, slong prec)
{
    //函数中所用变量
    arb_t s,a,t,w,v;
    
    arb_init(s);
    arb_init(a);
    arb_init(t);
    arb_init(w);
    arb_init(v);
    
    //前面部分
    //前系数
    arb_sqr(v,xi,prec);
    arb_sub_si(v,v,3,prec);
    arb_mul(v,v,xi,prec);
    arb_div_si(v,v,2,prec);
    
    //括号内
    arb_set_str(s,"2.5",prec); // 5/2
    arb_sqrt(s,s,prec);
    arb_mul(t,s,xi,prec);// t 后面多次用到 t=sqrt(5/2)*xi
    arb_div_si(s,t,2,prec);
    arb_hypgeom_erf(s,s,prec); //误差函数
    arb_hypgeom_erf(w,t,prec);
    arb_add(s,s,w,prec);
    
    arb_mul(v,v,s,prec);
    
    
    //后面括号前部分
    arb_set_str(s,"1.6",prec);// 8/5
    arb_sqr(w,xi,prec);
    arb_mul_si(w,w,31,prec);
    arb_div_si(w,w,4,prec);
    arb_add(s,s,w,prec);
    
    arb_sqr(w,t,prec);
    arb_div_si(w,w,4,prec);
    arb_neg(w,w);
    arb_exp(w,w,prec);
    
    arb_mul(s,s,w,prec);
    
    //后面括号后部分
    arb_sqr(w,t,prec);
    arb_div_si(w,w,5,prec);
    arb_set_str(a,"1.6",prec); // 8/5
    arb_sub(w,w,a,prec);
    
    arb_sqr(a,t,prec);
    arb_neg(a,a);
    arb_exp(a,a,prec);
    
    arb_mul(w,w,a,prec);
    
    arb_add(s,s,w,prec); //后面括号中的两部分加起来
    
    //后面系数部分
    arb_mul_si(w,Pi,5,prec);
    arb_inv(w,w,prec);
    arb_mul_si(w,w,2,prec);
    arb_sqrt(w,w,prec);
    
    arb_mul(s,s,w,prec);
    
    //两部分相加
    arb_add(res,v,s,prec);
    
    arb_clear(s);
    arb_clear(a);
    arb_clear(t);
    arb_clear(w);
    arb_clear(v);
    
    return 0;
}

//峰数密度相关函数 (P_1)^3(ν,ξ)
int N_pk_P_1_3(arb_t res, const arb_t nu, const arb_t xi,slong prec)
{
    //函数中所用变量
    arb_t s,t,w;
    
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //前面系数部分
    arb_sqr(s,Gamma_3,prec); //1-（Gamma_3）^2
    arb_neg(s,s);
    arb_add_si(s,s,1,prec);
    
    arb_sqrt(t,s,prec);
    arb_mul(t,t,Pi_2,prec);
    arb_inv(t,t,prec);
    
    //后面指数部分
    arb_mul(w,nu,Gamma_1,prec);
    arb_sub(w,xi,w,prec);
    arb_sqr(w,w,prec);
    arb_div(w,w,s,prec); // s 用掉
    
    arb_sqr(s,nu,prec);
    arb_add(s,s,w,prec);
    arb_div_si(s,s,2,prec);
    arb_neg(s,s);
    arb_exp(s,s,2*prec);
    
    arb_mul(res,t,s,prec);
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    return 0;
}


//峰的数密度，依照高斯类型来计算的 n_pk(mu_2,k_3)
//一般情问下，是mu_2和k_3的函数，需通过积分把k_3积掉
//在delta谱的情况下，大为化简，由于P_1^3 里出现delta函数，可以得到 n_pk(mu_2)的解析式

int Peak_number_density(arb_t res, const arb_t mu, const arb_t k, slong prec)
{
    // mu -> µ_2   k -> κ_3
    //函数中所用变量
    arb_t s,t,w,v;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(v);
    
    //跟功率谱有关，对于delta型式的功率谱，可解析求解
    //功率谱类型判断
    switch(Power_spectrum_type) 
    {
        case lognormal_type :
            
            //前面系数部分
            arb_one(s);
            arb_mul_si(s,s,3,prec); // 3
            arb_div_si(t,s,2,prec); // 3/2
            
            arb_pow(w,s,t,prec);//系数分子
            arb_mul_si(w,w,2,prec);
            arb_pow(s,Pi_2,t,prec);//系数分母
            arb_div(s,w,s,prec);
            
            //系数中 σ 部分
            arb_mul(s,s,Sigma_2_square,prec);
            arb_sqrt(v,Sigma_4_square,prec);
            arb_mul(s,s,Sigma_4_square,prec);
            arb_mul(s,s,v,prec);
            
            arb_div(s,s,Sigma_1_square,prec);
            arb_div(s,s,Sigma_1_square,prec);
            arb_sqrt(v,Sigma_3_square,prec);
            arb_div(s,s,Sigma_3_square,prec);
            arb_div(s,s,v,prec); //系数  σ 部分 完
            
            arb_mul(s,s,mu,prec);
            arb_mul(s,s,k,prec);
            
            arb_sqrt(v,Sigma_2_square,prec); //第一个娈量
            arb_div(v,v,Sigma_1_square,prec);
            arb_mul(v,v,mu,prec);
            
            arb_sqr(w,k,prec);//第二个娈量
            arb_mul(w,w,v,prec);
            
            N_pk_P_1_3(t,v,w,prec);//后面函数部分
            N_pk_f_xi(v,w,prec);
            arb_mul(t,t,v,prec);
            
            arb_mul(res,s,t,prec);
            
            break;
            
        case delta_type :
            //在delta谱的情况下，大为化简，由于P_1^3 里出现delta函数，可以得到 n_pk(mu_2)的解析式
            //此时，已与变量 k_3 无关了
            //前面常系数
            arb_inv(s,Pi_2,prec);
            arb_mul_si(s,s,3,prec);
            arb_sqrt(t,s,prec);
            arb_mul(w,s,t,prec);
            
            //中间部分
            arb_pow_ui(s,K_star,3,prec);
            arb_mul(w,w,s,prec);
            
            arb_sqrt(s,Power_A,prec);
            arb_div(s,mu,s,prec);
            N_pk_f_xi(t,s,prec);
            arb_mul(w,w,t,prec);
            
            arb_mul(s,Pi_2,Power_A,prec);
            arb_sqrt(s,s,prec);
            arb_inv(s,s,prec);
            arb_mul(w,w,s,prec);
            
            //尾部指数
            arb_mul_si(s,Power_A,2,prec);
            arb_sqr(t,mu,prec);
            arb_div(s,t,s,prec);
            arb_neg(s,s);
            arb_exp(s,s,prec);
            
            arb_mul(res,w,s,prec);
            
            break;
            
        default :
            printf("Peak_Theory -> number_density -> Peak_number_density -> power_spectrum_type 有误\n");
            exit(1);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(v);
    
    return 0;
}



//计算质量为M的原初黑洞数密度
int PBH_number_density_M(arb_t res,const arb_t M, slong prec)
{
    //函数中所用变量
    arb_t s,w,t_mu;
    arb_init(s);
    arb_init(w);
    arb_init(t_mu);
    
    // dln(M)/dµ 与k_3 无关
    //首先求 M 对应的 µ 值
    Horizon_reentry_mu_M(t_mu, M, prec); //这一步会有 R_MAX 保存
    
    //arb_set_str(t_mu,"0.01",prec); //测试
    
    
    //跟功率谱有关，对于delta型式的功率谱，可解析求解
    //功率谱类型判断
    switch(Power_spectrum_type) 
    {
        case lognormal_type :
            //不能解析求解，需要积分
            
            //对 n_pk(mu_2,k_3) 的 k_3 进行积分
            int interior_PBH_number_density_M(arb_t res, const arb_t k, void * param, const slong order, slong prec)
            {
                // µ_2 已由参数传入
                arb_t mu;
                arb_init(mu);
                
                arb_set(mu, param);
                
                Peak_number_density(res,mu,k,prec);
                arb_clear(mu);
                return 0;
            }
            
            int ret_judge=0;
            ret_judge=Integration_arb(s, interior_PBH_number_density_M, t_mu, 0,
                                                Int_n_pk_k_3_min, Int_n_pk_k_3_max,Int_n_pk_k_3_precision,
                                                Integration_iterate_min,Integration_iterate_max, prec);
            if(ret_judge==1)
            {
                printf("PBH_number_density_M \t 达到最大迭代次数\n");
            }
            break;
            
        case delta_type :
            //可解析求解，不需要对 k_3 进行积分
            
            Peak_number_density(s,t_mu,t_mu,prec); // 后一个mu作为k的补充参数，这里无实际作用
            
            break;
        default:
            printf("Peak_Theory -> number_density -> PBH_number_density_M -> power_spectrum_type 有误\n");
            exit(1);
    }
    
    
    Horizon_reentry_D_ln_M_to_mu(w, t_mu, prec); //可利用上面保存的 R_MAX
    arb_inv(w,w,prec); //取倒数
    
    arb_mul(res,s,w,prec);
    
    arb_clear(s);
    arb_clear(w);
    arb_clear(t_mu);
    
    return 0;
}

