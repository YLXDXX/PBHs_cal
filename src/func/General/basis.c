#include "basis.h"
#include <stdlib.h>

//窗口函数
int Power_spectra_window_function_k(arb_t res, const arb_t k, const arb_t R,
                                    const enum WINDOW_FUNC_TYPE w_type, slong prec)
{
    arb_t s,t,k_R;
    arb_init(s);
    arb_init(t);
    arb_init(k_R);
    
    //注意，这里的k传入时，是原本的k，没有取对数，传参时需注意
    
    //功率谱类型判断
    switch(w_type) 
    {
        case Real_space_top_hat : //形状参考 1905.01230
            //位置空间 W(r,R)=3/(4*π*R^3)*Θ(R-r)
            //傅氏空间 W(k,r)=3*( sin(k*R)-k*R*cos(k*R)/(k*R)^3 )
            
            arb_mul(k_R,k,R,prec); //k*R
            
            arb_cos(s,k_R,prec); //分子
            arb_mul(s,s,k_R,prec);
            arb_sin(t,k_R,prec);
            arb_sub(s,t,s,prec);
            arb_mul_ui(s,s,3,prec);
            
            arb_pow_ui(t,k_R,3,prec); //分母
            
            arb_div(res,s,t,prec);
            
            break;
        case Fourier_space_top_hat :
            //位置空间 W(r,R)= 1905.01230 (25)
            //傅氏空间 W(k,r)=Θ( 2.744/R - k )
            
            arb_set_str(s,"2.744",prec);
            arb_div(s,s,R,prec);
            arb_sub(s,s,k,prec);
            
            Heaviside_Theta_function(t,s,prec);
            
            arb_set(res,t);
            
            break;
        case Gaussian_hat :
            //位置空间 W(r,R)=1/(π*R^2)^(3/2)*Exp(-r^2/R^2)
            //傅氏空间 W(k,r)=Exp(- (k*R)^2/4 )
            
            arb_mul(k_R,k,R,prec); //k*R
            arb_sqr(s,k_R,prec);
            arb_div_ui(s,s,4,prec);
            arb_neg(s,s);
            
            arb_exp(res,s,prec);
            
            break;
        case Natural_hat_xx : //这后面几种，是为使用compaction function 计算准备的
                            //这里，形式上与前面对齐，求方差都为
                            //σ = ∫d(ln k) * W^2(k,R) * T^2(k,η) * P(k)
            //位置空间 W(r,R)=??
            //傅氏空间 W(k,r)=(k*r) * dj_0/d_z(kr), j_0(z)=sin(z)/z
            //dj_0/d_z=[z*cos(z)-sin(z)]/z^2
            
            arb_mul(k_R,k,R,prec); //k*R
            
            arb_cos(s,k_R,prec); //分母
            arb_mul(s,s,k_R,prec);
            arb_sin(t,k_R,prec);
            arb_sub(s,s,t,prec);
            
            arb_div(res,s,k_R,prec); //与前面约了一个k_R
            
            break;
        case Natural_hat_yy :
            //位置空间 W(r,R)=1/(4*π*R^2)*δ(r-R), δ(x) 是 delta 函数
            //傅氏空间 W(k,r)=j_0(k*R)
            
            arb_mul(k_R,k,R,prec); //k*R
            
            arb_sin(s,k_R,prec);
            arb_div(res,s,k_R,prec);
            
            break;
        default :
            printf("General -> basis -> 中 Power_spectra_window_function_k 类型有误\n");
            exit(1);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_init(k_R);
    
    return 0;
}



//线性转移函数
int Power_spectra_linear_transfer_function(arb_t res, const arb_t k, const arb_t eta, slong prec)
{
    arb_t x,s,t;
    
    arb_init(x);
    arb_init(s);
    arb_init(t);
    
    //注意，这里的k传入时，是原本的k，没有取对数，传参时需注意
    
    arb_mul(x,k,eta,prec); //x=kη/√3
    arb_sqrt_ui(s,3,prec);
    arb_div(x,x,s,prec);
    
    //T(x)=3*(sin(x)-x*cos(x))/x^3
    
    arb_cos(s,x,prec); //sin(x)-x*cos(x)
    arb_mul(s,s,x,prec);
    arb_sin(t,x,prec);
    arb_sub(s,t,s,prec);
    
    arb_pow_ui(t,x,3,prec);//x^3
    arb_div(s,s,t,prec);
    
    arb_mul_ui(res,s,3,prec);
    
    arb_clear(x);
    arb_clear(s);
    arb_clear(t);
    
    return 0;
}


//功率谱辅助函数
//这里，不同k对应的势能，关系后期推导，先直接用 V_0 近似
//Help_V_0b函数暂不使用
static void Help_V_0(arb_t res, const arb_t k, slong prec)
{
    //这里传入的k没有取对数
    
    arb_t s,t,w,q;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(q);
    
    //2207.11910 (4.1)
    //第一部分
    arb_mul_ui(t,Upward_step_spectra_epsilon_S,2,prec);
    arb_sqrt(t,t,prec);
    
    arb_sub(s,Upward_step_spectra_k_s,k,prec);
    Heaviside_Theta_function(w,s,prec);
    arb_mul(t,t,s,prec);
    arb_mul(t,t,w,prec);
    arb_add_ui(t,t,1,prec);
    
    arb_sub(s,Upward_step_spectra_k_c,k,prec);
    Heaviside_Theta_function(w,s,prec);
    arb_mul(t,t,Upward_step_spectra_V_0,prec);
    arb_mul(t,t,w,prec);
    
    //第二部分
    arb_mul_ui(s,Upward_step_spectra_epsilon_V,2,prec);
    arb_sqrt(s,s,prec);
    arb_sub(q,Upward_step_spectra_k_c,k,prec);
    arb_mul(s,s,q,prec);
    arb_add_ui(s,s,1,prec);
    
    arb_sqr(w,q,prec);
    arb_mul(w,w,Upward_step_spectra_eta_V,prec);
    arb_div_ui(w,w,2,prec);
    arb_add(s,s,w,prec);
    
    arb_neg(q,q);
    Heaviside_Theta_function(w,q,prec);
    arb_mul(s,s,w,prec);
    
    arb_add(w,Upward_step_spectra_V_0,Upward_step_spectra_Delta_V,prec);
    arb_mul(s,s,w,prec);
    
    arb_add(res,t,s,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(q);
}

static void Help_H_2_pi_2(arb_t res, const arb_t k, slong prec) //(H/2π)^2
{
    //这里传入的k没有取对数
    
    //[H/(2π*m_pl)]^2=(3*H^2 * m_pl^2)/(12*π^2*m_pl^4)
    //3H^2 * m_pl^2 ≈ V_0
    //由于在 V_0 的值中带了 m_pl^4 这个参数，可以约去
    //所以，在后续的计算中设 m_pl=1
    
    arb_t s,t,V_0;
    arb_init(s);
    arb_init(t);
    arb_init(V_0);
    
    //[H/(2π)]^2=V_0/(12*π^2)
    if ( 0 )
    {
        //这里，不同k对应的势能，关系后期推导，先直接用 V_0 近似 
        Help_V_0(V_0,k,prec);//不同的k对应于不同的势能
    }else
    {
        arb_set(V_0,Upward_step_spectra_V_0);
    }
    
    arb_sqr(t,Pi,prec);
    arb_mul_ui(t,t,12,prec);
    
    arb_div(res,V_0,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(V_0);
}

//由于 k_s 左右，有 k^4 关系
//通过拟合 10^-2*k_s 到 k_c 间 k^4 关系，给定 k_c 后，求 k_s
//此处， k^4 关系处理，不够完善，后期处理
void Upward_step_power_spectrum_k_c_to_k_s(arb_t k_s, const arb_t k_c, slong prec)
{
    arb_t P_kc,P_2ks,s,t,w;
    arb_init(P_kc);
    arb_init(P_2ks);
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //求P(10-2k_s)
    arb_div_ui(t,k_c,1E6,prec);
    Help_H_2_pi_2(s,t,prec); //传入未取对数的 k 
    
    arb_div(s,s,Upward_step_spectra_epsilon_S,prec);
    arb_div_ui(P_2ks,s,2,prec);
    arb_log(P_2ks,P_2ks,prec); //取对数
    
    //求P(k_c)
    Help_H_2_pi_2(s,k_c,prec); //传入未取对数的 k
    
    arb_sqr(w,Upward_step_spectra_g,prec);
    arb_mul(t,w,Upward_step_spectra_h,prec);
    arb_div_ui(t,t,6,prec);
    arb_neg(t,t);
    arb_add_ui(t,t,1,prec);
    arb_sqr(t,t,prec);
    arb_mul(s,s,t,prec);
    
    arb_mul(t,Upward_step_spectra_epsilon_V,w,prec);
    arb_mul_ui(t,t,2,prec);
    
    arb_div(P_kc,s,t,prec);
    arb_log(P_kc,P_kc,prec); //取对数
    
    //k^4变化关系
    //取对数后，为线性关系
    
    arb_log(t,k_c,prec);
    
    arb_sub(s,P_kc,P_2ks,prec);
    arb_div_ui(s,s,4,prec);
    
    arb_sub(w,t,s,prec);
    arb_exp(s,w,2*prec);
    arb_mul_ui(k_s,s,1E2,prec);
    
    arb_clear(P_kc);
    arb_clear(P_2ks);
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);    
}


//功率谱 
int power_spectrum(arb_t res, const arb_t k, slong prec)
{
    //函数中所用变量
    arb_t s,t,w;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //功率谱类型判断
    switch(Power_spectrum_type) 
    {
        case lognormal_type :
            
            //log-normal形式
            // P(k)=A/sqrt(2*Pi*sigma^2)*exp( -( ln(k)-ln(k_star) )^2/(2*sigma^2) )
            
            //arb_log(t, k, prec); //以 k 作为变量
            //arb_sub(t,t,Ln_K_star,prec); // 以 k 或 ln(k) 这里不一样
            
            arb_sub(t,k,Ln_K_star,prec); //以 ln(k) 作为变量
            
            //指数部分
            arb_sqr(t,t,prec);
            arb_sqr(s,Power_sigma,prec);
            arb_mul_si(s,s,2,prec);
            arb_div(t,t,s,prec);
            arb_neg(t,t);
            arb_exp(t,t,prec);
            
            //系数部分
            arb_mul(t,t,Power_A,prec); // A/sqrt(2*Pi*sigma^2)=A/[sqrt(2*Pi)*Δ]
            arb_sqrt(s,Pi_2,prec);
            arb_mul(s,s,Power_sigma,prec);
            arb_div(res,t,s,prec);
            
            
            break;
            
        case power_law_type : //P(k)=2*π^2*A*k^{-3} --> P(ln(k))=2*π^2*A*[e^(ln(k))]^-3
            //这里以ln(k)作为变量，传入的即是ln(k)的值
            // P(k)=2*π^2*A*e^(-3k)
            arb_sqr(t,Pi,prec);
            arb_mul_ui(t,t,2,prec);
            arb_mul(t,t,Power_A,prec);
            
            arb_exp(s,k,prec); //k^{-3}=[e^ln(k)]^(-3)
            arb_pow_ui(s,s,3,prec);
            arb_inv(s,s,prec);
            
            arb_mul(res,t,s,prec);
            
            break;
            
        case box_type :
            //注意，这里是以ln(k)作为自变量的，传进来的值就是ln(k)的值
            //Box 区间[s,t],区间内为A,区间外为0
            arb_div_ui(w,Box_width,2,prec);
            
            arb_sub(s,Box_K_middle,w,prec); //s=K_middle - width/2
            arb_log(s,s,prec); //ln(k)作变量，取对数
            
            arb_add(t,Box_K_middle,w,prec); //t=K_middle + width/2
            arb_log(t,t,prec); //ln(k)作变量，取对数
            
            if( arb_gt(k,s) && arb_lt(k,t) ) // k>s && k<t
            {
                arb_set(res,Power_A);
            }else
            {
                arb_zero(res);
            }
            
            break;
            
        case broken_power_law_type :
            //broken power law 功率谱
            //P(k)=A*(α+β)^γ/[ β*(k/k_star)^(-α/γ) + α*(k/k_star)^(β/γ) ]^γ
            //注意，这里是以ln(k)作为自变量，传入的是ln(k)
            
            arb_exp(w,k,prec); //k=e^ln(k)
            arb_div(w,w,K_star,prec); // k/k_star
            
            
            //分母部分 [ β*(k/k_star)^(-α/γ) + α*(k/k_star)^(β/γ) ]^γ
            arb_div(t,BPL_alpha,BPL_gamma,prec); //β*(k/k_star)^(-α/γ)
            arb_neg(t,t);
            arb_pow(s,w,t,prec); 
            arb_mul(s,s,BPL_beta,prec);
            
            arb_div(t,BPL_beta,BPL_gamma,prec); //α*(k/k_star)^(β/γ)
            arb_pow(t,w,t,prec); 
            arb_mul(t,t,BPL_alpha,prec);
            
            arb_add(t,t,s,prec);
            arb_pow(t,t,BPL_gamma,prec); 
            
            //分子部分A*(α+β)^γ
            arb_add(s,BPL_alpha,BPL_beta,prec);
            arb_pow(s,s,BPL_gamma,prec);
            arb_mul(s,s,Power_A,prec);
            
            arb_div(res,s,t,prec);
            
            break;
            
        case link_cmb_type :
            //连接CMB功率谱，一个三段函数
            
            arb_exp(w,k,prec);  //传入的是k=e^ln(k)
            
            if( arb_lt(w,Link_CMB_K_t) ) //k<k_t
            {
                //第一段 A_s*(k/k_star)^(n_s-1)
                arb_div(s,w,Link_CMB_K_star,prec);
                arb_sub_ui(t,Link_CMB_n_s,1,prec);
                arb_pow(s,s,t,prec);
                
                arb_mul(res,s,Link_CMB_A_s,prec);
                
            }else if( arb_lt(w,Link_CMB_K_m) ) //k<k_m
            {
                //第二段 A_b*(k/k_star)^(n_b-1)
                arb_div(s,w,Link_CMB_K_star,prec);
                arb_sub_ui(t,Link_CMB_n_b,1,prec);
                arb_pow(s,s,t,prec);
                
                arb_mul(res,s,Link_CMB_A_b,prec);
                
            }else
            {
                //第三段 P_m(k)
                arb_set(res,Link_CMB_P_m);
            }
            
            break;
            
        case upward_step_spectra_type :
            //注意，这里是以ln(k)作为自变量的，传进来的值就是ln(k)的值
            //2207.11910 (4.13)
            
            arb_set_str(t,"1E-2",prec); //ln(1E-2*k_s)
            arb_log(t,t,prec);
            arb_add(t,t,Upward_step_spectra_ln_k_s,prec);
            
            if ( arb_lt(k,t) ) //原式为 k ≪ k_s，这里使用 k < 1E-2*k_s
            {
                arb_exp(t,k,prec);
                Help_H_2_pi_2(s,t,prec); //传入未取对数的 k 
                
                arb_div(s,s,Upward_step_spectra_epsilon_S,prec);
                arb_div_ui(res,s,2,prec);
                
            }else if ( arb_le(k,Upward_step_spectra_ln_k_c) ) // k < k_c
            {
                //在 1E-2*k_s 到 k_c 这一段，有一个 k^4 变化上升的规律
                
                arb_exp(t,Upward_step_spectra_ln_k_c,prec);
                Help_H_2_pi_2(s,t,prec); //传入未取对数的 k
                
                arb_sqr(w,Upward_step_spectra_g,prec);
                arb_mul(t,w,Upward_step_spectra_h,prec);
                arb_div_ui(t,t,6,prec);
                arb_neg(t,t);
                arb_add_ui(t,t,1,prec);
                arb_sqr(t,t,prec);
                arb_mul(s,s,t,prec);
                
                arb_mul(t,Upward_step_spectra_epsilon_V,w,prec);
                arb_mul_ui(t,t,2,prec);
                
                arb_div(s,s,t,prec); //这里求得 P(k_c)
                arb_log(s,s,prec); //取对数
                
                //P(k)=P(k_c)-(ln(k_c)-ln(k))*4
                arb_sub(t,Upward_step_spectra_ln_k_c,k,prec);
                arb_mul_ui(t,t,4,prec);
                arb_sub(s,s,t,prec);
                arb_exp(res,s,prec); //对数结果，取指数还原
                
            }else // k > k_c
            {
                arb_pow_ui(w,Upward_step_spectra_g,4,prec);
                arb_neg(t,w);
                arb_add_ui(t,t,1,prec);
                
                arb_exp(s,k,prec);
                arb_mul(s,s,Upward_step_spectra_tau_c,prec);
                arb_mul_ui(s,s,2,prec);
                arb_cos(s,s,prec);
                
                arb_mul(t,t,s,prec);
                arb_add_ui(t,t,1,prec);
                arb_add(t,t,w,prec);
                
                arb_sqr(w,Upward_step_spectra_g,prec);
                arb_mul_ui(w,w,4,prec);
                arb_mul(w,w,Upward_step_spectra_epsilon_V,prec);
                
                arb_div(w,t,w,prec);
                
                
                
                arb_exp(t,k,prec);
                Help_H_2_pi_2(s,t,prec); //传入未取对数的 k
                
                arb_mul(res,w,s,prec);
            }
            
            break;
            
        case delta_type :
            
            //delta形式
            //P(k)=A*delta( ln(k)-ln(k_star) )
            //可用 log-normal 形式取极限 σ → 0 得到
            
            arb_t tem_sigma;
            arb_init(tem_sigma);
            arb_set_str(tem_sigma,"1E-13",100);
            
            arb_sub(t,k,Ln_K_star,prec); //以 ln(k) 作为变量
            
            //指数部分
            arb_sqr(t,t,prec);
            arb_sqr(s,tem_sigma,prec);
            arb_mul_si(s,s,2,prec);
            arb_div(t,t,s,prec);
            arb_neg(t,t);
            arb_exp(t,t,prec);
            
            //系数部分
            arb_mul(t,t,Power_A,prec);
            arb_sqrt(s,Pi_2,prec);
            arb_mul(s,s,Power_sigma,prec);
            arb_div(res,t,s,prec);
            
            arb_clear(tem_sigma);
            
            break;
            
        default :
            printf("General -> basis -> power_spectrum 中 power_spectrum_type 有误\n");
            exit(1);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    
    return 0;
}




