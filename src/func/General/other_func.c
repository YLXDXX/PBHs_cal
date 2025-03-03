#include "other_func.h"
#include <stdlib.h>

//将compaction func C 的值转为 C_l 的值
int Trans_C_to_C_l(arb_t res, const arb_t x, slong prec)
{
    arb_t s,t;
    
    arb_init(s);
    arb_init(t);
    
    //二次方程，有两个根，考虑到第一型扰动，取其中的一个根
    //C_l_th=4/3*( 1-sqrt(1-3/2*x) )
    
    arb_one(t); 
    arb_mul_ui(s,t,3,prec); //sqrt(1-3/2*x)
    arb_div_ui(s,s,2,prec);
    arb_mul(s,s,x,prec);
    arb_neg(s,s);
    arb_add_ui(s,s,1,prec);
    arb_sqrt(s,s,prec);
    
    arb_neg(s,s);
    arb_add_ui(s,s,1,prec);
    
    arb_mul_ui(s,s,4,prec);
    arb_div_ui(res,s,3,prec);
    
    arb_clear(s);
    arb_clear(t);
    
    return 0; 
}


static int interior_Get_PK_mu_max(arb_t res, const arb_t mu, void * zeta_k, const slong order, slong prec)
{
    arb_t s;
    arb_init(s);
    
    //与 μ 和 k 均有关
    arb_set(PT_mu,mu);
    arb_set(PT_k,zeta_k);
    
    //通过利用 C_l 的最大值 4/3 来得到 μ 的最大值
    //C_l=-4/3*r*ζ' --> r*ζ'的最大值为 -1
    
    zeta_profile_n(s, R_MAX, 1, prec);
    arb_mul(s,s,R_MAX,prec);
    
    arb_add_ui(res,s,1,prec);
    
    arb_clear(s);
    
    return 0; 
}

//求出参数 μ 的上限
int Get_PK_mu_max(arb_t res, const arb_t zeta_k, slong prec)
{
    arb_t s,a,b,k;
    
    arb_init(s);
    arb_init(a);
    arb_init(b);
    arb_init(k);
    
    //当不考虑profile的简化时，与k参数相关
    arb_set(k,zeta_k);
    
    //找根的最小值，从 PT_mu_th 开始
    arb_set(a,PT_mu_th); 
    arb_set(b,Int_mu_max);
    arb_mul_ui(b,b,2,prec);
    
    Find_interval_root(s, interior_Get_PK_mu_max, k, 0,
                       a, b, Int_mu_precision,
                       2*Root_mu_num, Root_Normal, prec);
    arb_set(res,s);
    
    arb_clear(s);
    arb_clear(a);
    arb_clear(b);
    arb_clear(k);
    
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
            //参见 2109.00791 (A.15) 第一行
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
            //参见 2109.00791 (A.15) 第二行
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
        case beta_f_general_III :
            //自己文章中采用的转换公式
            //即将 2109.00791 (A.15) 第一行：f(m) = Ω_m/Ω_DM * T_{r_m}/T_eq * β(M) 中的
            //T_{r_m}/T_eq 利用波数和温度的关系替换「见 1812.00674 的 (A5)」
            //f(M)=Ω_m/Ω_DM * (g_{*,eq}/g_*)^(1/6) * (k_*/k_eq) * 1/[2*(√2-1)] * β(M)
            
            arb_div(t,Omega_M,Omega_DM,prec);
            arb_div(s,effective_g_star_eq,effective_g_star,prec);
            arb_one(w);
            arb_div_ui(w,w,6,prec);
            arb_pow(s,s,w,prec);
            arb_mul(t,t,s,prec);
            
            arb_div(s,K_star,K_scale_eq,prec); //k_*/k_eq
            arb_mul(t,t,s,prec);
            
            arb_one(s); // 1/[2*(√2-1)]
            arb_mul_ui(s,s,2,prec);
            arb_sqrt(s,s,prec);
            arb_sub_ui(s,s,1,prec);
            arb_mul_ui(s,s,2,prec);
            
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

//考虑所有 k 模式时，特征模式的求解
void Get_all_k_over_k_ch(arb_t k_ch_times_r_m, arb_t k_ch, const arb_t x_m, slong prec)
{
    arb_t t,s;
    arb_init(t);
    arb_init(s);
    
    //有两种方法，一种是利用δ谱近似来求，一种是利用进入视界的条件来求
    //这里的变量 x_m 仅在利用δ谱求时使用
    if(Get_k_ch_type==delta_approximation)
    {
        //利用δ近似
        //k_ch=x_m/r_m 
        arb_div(k_ch,x_m,R_MAX,prec);
        
    }else if (Get_k_ch_type==horizon_re_enter)
    {
        //利用进入视界的条件
        //k_ch=π/[r_m*e^(ζ_m)]
        zeta_profile_n(s,R_MAX,0,prec);
        arb_exp(s,s,prec);
        arb_mul(s,s,R_MAX,prec);
        
        arb_div(s, Pi, s, prec);
        if(Power_spectrum_type!=delta_type)
        {
            arb_mul(s,s,Delta_spectrum_reenter_coefficient_C,prec);
        }
        arb_set(k_ch,s);
    }else
    {
        //出错
        printf("连续谱考虑所有k模式，求特征模式方法 Get_k_ch_type 错误\n");
        exit(1);
    }
    
    arb_mul(k_ch_times_r_m,k_ch,R_MAX,prec);
    
    arb_clear(t);
    arb_clear(s);
}


