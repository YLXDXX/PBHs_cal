#include "basis.h"
#include <stdlib.h>

//在k空间的窗口函数
int window_function(arb_t res, const arb_t k, slong prec)
{
    arb_t s;
    arb_init(s);
    
    //功率谱类型判断
    switch(Power_spectrum_type) 
    {
        case lognormal_type :
            //高斯型的窗口函数，中心点在 k_star
            //exp[-(k-k_star)^2]
            
            //与功率谱等一同修改
            //arb_sub(s,k,K_star,prec); //以 k 为变量
            arb_exp(s,k,prec); //以 ln(k) 为变量， k=e^ln(k)
            arb_sub(s,s,K_star,prec);
            
            arb_sqr(s,s,prec);
            arb_neg(s,s);
            arb_exp(res,s,prec);
            
            break;
            
        default :
            printf("General -> basis -> window_function 中 power_spectrum_type 有误\n");
            exit(1);
    }
    
    arb_clear(s);
    
    return 0;
}



//线性转移函数
int Linear_transfer_function(arb_t res, const arb_t k, const arb_t eta, slong prec)
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
            
            //与窗口函数等一同修改
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
            arb_mul(t,t,Power_A,prec); // A/sqrt(2*Pi*sigma^2)=A/sqrt(2*Pi)*Δ
            arb_sqrt(s,Pi_2,prec);
            arb_mul(s,s,Power_sigma,prec);
            arb_div(res,t,s,prec);
            
            //加上窗口函数
            //window_function(s,k,prec);
            //arb_mul(res,t,s,prec);
            
            //加入转移函数 P_ζ(k,η)=T^2(k,η)*P_ζ(k)
            //用Peake theory来算阈值，加入转移函数的改变几乎可以忽略不计
            //在算方差的时候，再单独考虑转移函数的影响
            /*
            arb_set_str(t,"1.75878702981131e-13",prec);
            arb_exp(w,k,prec);
            Linear_transfer_function(s, w, t, prec);//转移函数需传入原本的未取对数的k
            arb_sqr(s,s,prec);
            arb_mul(res,res,s,prec);
            */
            
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
            arb_div(res,t,s,prec); //不需要窗口函数了
            
            //加上窗口函数
            //window_function(s,k,prec);
            //arb_mul(res,t,s,prec);
            
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




