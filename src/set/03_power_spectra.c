#include "header/03_power_spectra.h"


void Set_power_spectra(char* comd_argv, slong prec) // comd_argv 为命令行传递参数
{
    //log-normal谱
    arb_set_str(Power_A, "4E-3", prec); //功率谱振幅
    
    //功率谱掌宽，例如 log-normal 的
    arb_set_str(Power_sigma,"0.1",prec); // P(k)=A/sqrt(2*Pi*sigma^2)*exp( -( ln(k)-ln(k_star) )^2/(2*sigma^2) )
    //arb_set_str(Power_sigma,argv[1],prec);//命令行中读取
    
    
    //lognormal 功率谱参考尺度，单位 Mpc^{-1}
    arb_set_str(K_star,"1.56E13",prec); //参考值 K_star=1.56E13
    //arb_set_str(K_star,argv[1],prec);//命令行中读取
    //arb_exp(K_star,K_star,prec);
    
    
    //功率谱参考尺度，单位 Mpc^{-1}, 取对数
    arb_log(Ln_K_star,K_star,prec);
    
    
    //在log-normal谱的处理中，多处需要功率谱参与的积分
    //如：ζ_G(r)辅助变量、协方差阵、δ谱计算连续谱
    //积分区域跟谱宽Δ和峰位置k_*有关
    //功率谱掌宽的倍数
    arb_mul_ui(Log_normal_mul_sigma,Power_sigma, 6, prec); //动态积分范围取6σ，后面用
    
    
    //Box功率谱 中心点
    arb_set_str(Box_K_middle,"5",prec);
    
    //Box功率谱 宽度
    arb_set_str(Box_width,"1",prec);
    
    
    //broken power law功率谱使用 α
    //α 影响的是P(k_star)左边的部分，k<k_star，α越大，下降的越厉害
    arb_set_str(BPL_alpha,"20",prec);
    
    
    //broken power law功率谱使用 β 
    //β 影响的是P(k_star)右边的部分，k>k_star，β越大，下降的越厉害
    arb_set_str(BPL_beta,"20",prec);
    
    //broken power law功率谱使用 γ
    //γ 影响的是P(k)左右两边的转化，γ越大，转化越圆润，P(k_star)周围下降的越慢
    arb_one(BPL_gamma); //一般设 γ=1
    
    //link cmb 功率谱用
    arb_set_str(Link_CMB_K_t,"10",prec);
    arb_set_str(Link_CMB_K_m,"3E3",prec);
    arb_set_str(Link_CMB_A_s,"2.101E-9",prec);
    arb_set_str(Link_CMB_A_b,"1.64E-13",prec);
    arb_set_str(Link_CMB_K_star,"0.05",prec);
    arb_set_str(Link_CMB_P_m,"1E-6",prec);
    arb_set_str(Link_CMB_n_s,"0.9649",prec);
    arb_set_str(Link_CMB_n_b,"2.75",prec);
    
}
 
