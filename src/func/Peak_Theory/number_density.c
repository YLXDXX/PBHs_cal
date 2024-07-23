#include "number_density.h"
#include <stdlib.h>
#include <arb_hypgeom.h> //误差函数用

//
//主要参考2109.00791编写
//



//峰数密度计算相关辅助函数 f(ξ)
int N_pk_help_f_xi(arb_t res, const arb_t xi, slong prec)
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
    arb_one(s); // 5/2
    arb_mul_ui(s,s,5,prec);
    arb_div_ui(s,s,2,prec);
    
    arb_sqrt(s,s,prec);
    arb_mul(t,s,xi,prec);// t 后面多次用到 t=sqrt(5/2)*xi
    arb_div_si(s,t,2,prec);
    arb_hypgeom_erf(s,s,prec); //误差函数
    arb_hypgeom_erf(w,t,prec);
    arb_add(s,s,w,prec);
    
    arb_mul(v,v,s,prec);
    
    
    //后面括号前部分
    arb_one(s); // 8/5
    arb_mul_ui(s,s,8,prec);
    arb_div_ui(s,s,5,prec);
    
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
    arb_one(a); // 8/5
    arb_mul_ui(a,a,8,prec);
    arb_div_ui(a,a,5,prec);
    arb_sub(w,w,a,prec);
    
    arb_sqr(a,t,prec);
    arb_neg(a,a);
    arb_exp(a,a,prec);
    
    arb_mul(w,w,a,prec);
    
    arb_add(s,s,w,prec); //后面括号中的两部分加起来
    
    //后面系数部分
    arb_mul_si(w,Pi,5,prec);
    arb_ui_div(w,2,w,prec);
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

//峰数密度计算相辅助函数 (P_1)^n(ν,ξ)
int N_pk_help_P_1_n(arb_t res, const arb_t nu, const arb_t xi, const slong n, slong prec)
{
    //函数中所用变量
    arb_t s,t,w,gamma;
    
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(gamma);
    
    if(n==1)
    {
        arb_set(gamma,Gamma_1);
    }else if(n==3)
    {
        arb_set(gamma,Gamma_3);
    }else
    {
        printf("N_pk_help_P_1_n 阶数 n 输入有误\n");
        exit(1);
    }
    
    //前面系数部分
    arb_sqr(s,gamma,prec); //1-（Gamma_3）^2
    arb_neg(s,s);
    arb_add_si(s,s,1,prec);
    
    arb_sqrt(t,s,prec);
    arb_mul(t,t,Pi_2,prec);
    
    //后面指数部分
    arb_mul(w,nu,Gamma_1,prec);//原文公式此处漏写了下标1
    arb_sub(w,xi,w,prec);
    arb_sqr(w,w,prec);
    arb_div(w,w,s,prec); // s 用掉 
    
    arb_sqr(s,nu,prec);
    arb_add(s,s,w,prec);
    arb_div_si(s,s,2,prec);
    arb_neg(s,s);
    arb_exp(s,s,2*prec);
    
    arb_div(res,s,t,prec);
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(gamma);
    return 0;
}


//峰的数密度，依照高斯类型来计算的 n_pk(mu,k)
//一般情问下，是mu和k的函数，需通过积分把k积掉
//在delta谱的情况下，大为化简，由于P_1^3 里出现delta函数，可以得到 n_pk(mu)的解析式

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
        case delta_type :
            //没取梯度参见 1906.06790 (25)，
            //取梯度参见： 2109.00791 (3.11)
            //σ_n^2=A*k_star^(2n) --> σ_0^2=A，两者的结果一样
            //在delta谱的情况下，大为化简，由于P_1^n 里出现delta函数，可以得到 n_pk(mu)的解析式
            //此时，已与变量 k 无关了
            
            //前面常系数
            arb_ui_div(s,3,Pi_2,prec);
            arb_sqrt(t,s,prec);
            arb_mul(w,t,s,prec);
            
            //中间部分
            arb_pow_ui(s,K_star,3,prec);
            arb_mul(w,w,s,prec);
            
            arb_sqrt(s,Power_A,prec);
            arb_div(s,mu,s,prec);
            N_pk_help_f_xi(t,s,prec);
            arb_mul(w,w,t,prec);
            
            arb_mul(s,Pi_2,Power_A,prec);
            arb_sqrt(s,s,prec);
            arb_div(w,w,s,prec);
            
            //尾部指数
            arb_mul_si(s,Power_A,2,prec);
            arb_sqr(t,mu,prec);
            arb_div(s,t,s,prec);
            arb_neg(s,s);
            arb_exp(s,s,prec);
            
            arb_mul(res,w,s,prec);
            
            break;
            
         default :
            
            int n;
            if(Peak_theory_source_zeta_gradient==true) //peak theory计算中，统计量是否取ζ的梯度
            {
                n=3; //取了梯度
                
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
                
                N_pk_help_P_1_n(t,v,w,n,prec);//后面函数部分
                N_pk_help_f_xi(v,w,prec);
                arb_mul(t,t,v,prec);
                
                arb_mul(res,s,t,prec);
            }else
            {
                n=1; //没取梯度
                
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
                arb_sqrt(v,Sigma_0_square,prec);
                arb_div(s,s,v,prec);
                
                arb_sqrt(v,Sigma_1_square,prec);
                arb_div(s,s,v,prec);
                arb_div(s,s,Sigma_1_square,prec); //系数  σ 部分 完
                
                arb_mul(s,s,mu,prec);
                arb_mul(s,s,k,prec);
                
                arb_sqrt(v,Sigma_0_square,prec); //第一个娈量
                arb_div(v,mu,v,prec);
                
                arb_sqr(w,k,prec);//第二个娈量
                arb_mul(w,w,mu,prec);
                arb_sqrt(t,Sigma_2_square,prec);
                arb_div(w,w,t,prec);
                
                N_pk_help_P_1_n(t,v,w,n,prec);//后面函数部分
                N_pk_help_f_xi(v,w,prec);
                arb_mul(t,t,v,prec);
                
                arb_mul(res,s,t,prec);
            }
            
            break;
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(v);
    
    return 0;
}

//对 n_pk(mu,k) 的 k 进行积分，简化版，此时μ与k无关
static int interior_PBH_number_density_M_profile_simplify(arb_t res, const arb_t k, void * para_mu, const slong order, slong prec)
{
    arb_t s,t,mu;
    
    arb_init(s);
    arb_init(t);
    arb_init(mu);
    
    arb_set(mu,para_mu); //mu与k无关，由参数传入
    
    Peak_number_density(s,mu,k,prec);
    
    arb_set(res,s);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(mu);
    
    return 0;
}

//对 n_pk(mu,k) 的 k 进行积分，完整版
static int interior_PBH_number_density_M(arb_t res, const arb_t k, void * M, const slong order, slong prec)
{
    // M 已由参数传入
    arb_t s,t,mu,mu_th,mu_max;
    
    arb_init(s);
    arb_init(t);
    arb_init(mu);
    arb_init(mu_th);
    arb_init(mu_max);
    
    
    //首先求 M 对应的 µ 值
    Horizon_reentry_M_to_mu(mu, M, k, prec); //这一步会有 R_MAX 保存
    
    if(PT_threshold_simplify==true)
    {
        //此时不用再重新求μ的threshold和上限
        arb_set(mu_th,PT_mu_th);
        arb_set(mu_max,PT_mu_max);
        
    }else
    {
        //对每个 k 都重新求一遍
        Find_PT_Mu_th(mu_th, k, prec); //不同的 k 对应于不同的 threshold 
        Get_PK_mu_max(mu_max, PT_k, prec); //不同的 k 对应于不同的上限 
    }
    
    
    if( arb_lt(mu,mu_th) || arb_gt(mu,mu_max) ) //μ取值范围外，直接设为零
    {
        //此处，定义域无法显式写出，可以通过这种方法来积
        //只要积分区间足够大，覆盖相应的定义域即可
        
        arb_zero(res);
        
        arb_clear(s);
        arb_clear(t);
        arb_clear(mu);
        arb_clear(mu_th);
        arb_clear(mu_max);
        
        return 0;
    }else
    {
        Peak_number_density(s,mu,k,prec);
        Horizon_reentry_derivative_ln_M_mu(t, mu, k, prec); //可利用上面保存的 R_MAX
        arb_inv(t,t,prec); //取倒数
        
        arb_mul(res,s,t,prec);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(mu);
    arb_clear(mu_th);
    arb_clear(mu_max);
    
    return 0;
}


//计算质量为M的原初黑洞数密度
int PBH_number_density_M(arb_t res,const arb_t M, slong prec)
{
    //函数中所用变量
    arb_t s,w,t_mu,zeta_k;
    arb_init(s);
    arb_init(w);
    arb_init(t_mu);
    arb_init(zeta_k);
    
    int ret_judge=0;
    
    //PT_threshold_simplify=true
    //PT_profile_simplify==true
    
    if(PT_profile_simplify==true || Power_spectrum_type==delta_type) //包含 δ 谱
    {
        //上面判断两都的共同点是，μ的取值与k无关
        // dln(M)/dµ 也与 k 无关，此时计算较为简单
        
        //此时，k_3 为一常数，取其为对应的平均值
        //而 k 的取值，又要分有没有用 ζ 的梯度
        
        if(Peak_theory_source_zeta_gradient==true) //这里对于δ一样，但对于PT_profile_simplify不一样
        {
            //(k_3)^2=γ_3
            arb_sqrt(zeta_k,Gamma_3,prec);
        }else
        {
            //k_1=σ_1/σ_0
            arb_div(zeta_k,Sigma_1_square,Sigma_0_square,prec);
            arb_sqrt(zeta_k,zeta_k,prec);
        }
        
        
        //首先求 M 对应的 µ 值
        Horizon_reentry_M_to_mu(t_mu, M, zeta_k, prec); //这一步会有 R_MAX 保存
        
        if( arb_lt(t_mu,PT_mu_th) || arb_gt(t_mu,PT_mu_max) ) //μ取值范围外，直接设为零
        {
            arb_zero(res);
            
            arb_clear(s);
            arb_clear(w);
            arb_clear(t_mu);
            arb_clear(zeta_k);
            
            return 0;
        }else
        {
            if(Power_spectrum_type==delta_type) //δ的情况最简单，不用积分，直接代入计算
            {
                //积分化为常量积分
                Peak_number_density(s, t_mu, zeta_k, prec);
                
                Horizon_reentry_derivative_ln_M_mu(w, t_mu, zeta_k, prec); //可利用上面保存的 R_MAX
                arb_inv(w,w,prec); //取倒数
                arb_mul(res,s,w,prec);
                
            }else
            {
                //对于 PT_profile_simplify 的情况，还需要对于μ进行积分
                //但是，由于这里μ与k无关，可以大大加快积分的计算速度
                
                ret_judge=Integration_arb(s, interior_PBH_number_density_M_profile_simplify, t_mu, 0, //积分计算传递质量 M 
                                          Int_n_pk_k_min, Int_n_pk_k_max,Int_n_pk_k_precision,
                                          Integration_iterate_min,Integration_iterate_max, prec);
                
                Horizon_reentry_derivative_ln_M_mu(w, t_mu, zeta_k, prec); //可利用上面保存的 R_MAX
                arb_inv(w,w,prec); //取倒数
                
                arb_mul(res,s,w,prec);
            }
        }
        
    }else
    {
        //一般不能解析求解，计算较为复杂
        arb_set(w,M);
        
        ret_judge=Integration_arb(s, interior_PBH_number_density_M, w, 0, //积分计算传递质量 M 
                                  Int_n_pk_k_min, Int_n_pk_k_max,Int_n_pk_k_precision,
                                  Integration_iterate_min,Integration_iterate_max, prec);
        arb_set(res,s);
    }
    
    
    if(ret_judge==1)
    {
        printf("PBH_number_density_M \t 达到最大迭代次数\n");
    }
    
    arb_clear(s);
    arb_clear(w);
    arb_clear(t_mu);
    arb_clear(zeta_k);
    
    return 0;
}

