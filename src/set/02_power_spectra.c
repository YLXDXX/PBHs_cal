#include "header/02_power_spectra.h"
#include "../func/General/basis.h"

void Set_power_spectra(char* comd_argv, slong prec) // comd_argv 为命令行传递参数
{
    //log-normal谱
    arb_set_str(Power_A, "4E-3", prec); //功率谱振幅
    //arb_set_str(Power_A, "0.302", prec); //PTAs对log-normal的拟合
    //arb_set_str(Power_A, "0.211", prec); //PTAs对BPL的拟合
    //arb_set_str(Power_A, "0.164", prec); // upward step 拟合
    //arb_set_str(Power_A, comd_argv, prec);
    
    //功率谱掌宽，例如 log-normal 的
    arb_set_str(Power_sigma,"0.1",prec); // P(k)=A/sqrt(2*Pi*sigma^2)*exp( -( ln(k)-ln(k_star) )^2/(2*sigma^2) )
    //arb_set_str(Power_sigma,argv[1],prec);//命令行中读取
    //对于非δ、log-normal、BPL等中心对称谱，可选择某个合适的值
    //如，Peak thoery中会用到
    
    //lognormal 功率谱参考尺度，单位 Mpc^{-1}
    arb_set_str(K_star,"1.56E13",prec); //参考值 K_star=1.56E13
    //arb_set_str(K_star,"151356124.8",prec); //10^8.18 PTAs对log-normal的拟合
    //arb_set_str(K_star,"114815362.1",prec); //10^8.06 PTAs对BPL的拟合
    //arb_set_str(K_star,"6E7",prec); // upward step 拟合
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
    //arb_set_str(BPL_alpha,"4.52",prec); //PTAs对BPL的拟合
    arb_set_str(BPL_alpha,"4",prec); // upward step 拟合
    
    //broken power law功率谱使用 β 
    //β 影响的是P(k_star)右边的部分，k>k_star，β越大，下降的越厉害
    //arb_set_str(BPL_beta,"5.17",prec); //PTAs对BPL的拟合
    arb_set_str(BPL_beta,"6",prec); // upward step 拟合
    
    //broken power law功率谱使用 γ
    //γ 影响的是P(k)左右两边的转化，γ越大，转化越圆润，P(k_star)周围下降的越慢
    //arb_one(BPL_gamma); //一般设 γ=1
    //arb_set_str(BPL_gamma,"5.03",prec); //PTAs对BPL的拟合
    arb_set_str(BPL_gamma,"2",prec); // upward step 拟合
    
    //link cmb 功率谱用
    arb_set_str(Link_CMB_K_t,"10",prec);
    arb_set_str(Link_CMB_K_m,"3E3",prec);
    arb_set_str(Link_CMB_A_s,"2.101E-9",prec);
    arb_set_str(Link_CMB_A_b,"1.64E-13",prec);
    arb_set_str(Link_CMB_K_star,"0.05",prec);
    arb_set_str(Link_CMB_P_m,"1E-6",prec);
    arb_set_str(Link_CMB_n_s,"0.9649",prec);
    arb_set_str(Link_CMB_n_b,"2.75",prec);
    
    //upward step模型的原生功率谱, SR-USR-SR
    arb_set_str(Upward_step_spectra_epsilon_S,"2.551E-3",prec);
    arb_set_str(Upward_step_spectra_epsilon_V,"3.673E-7",prec);
    arb_set_str(Upward_step_spectra_eta_V,"7.143E-3",prec);
    arb_set_str(Upward_step_spectra_g,"0.007",prec);
    arb_set_str(Upward_step_spectra_h,"-5.5",prec);
    arb_set_str( Upward_step_spectra_V_0,"7E-10",prec); //7E-10*M_pl^4
    arb_set_str(Upward_step_spectra_Delta_V,"6.4E-16",prec); //6.4E-16*M_pl^4
    arb_set_str(Upward_step_spectra_Delta_phi_USR,"0.023",prec); //0.023*M_pl
    //arb_set_str(Upward_step_spectra_phi_s,"",prec);
    //arb_set_str(Upward_step_spectra_phi_c,"",prec);
    arb_set_str(Upward_step_spectra_k_c,"6E7",prec); //Mpc^-1
    //arb_set_str(Upward_step_spectra_k_s,"1E8",prec); //Mpc^-1
    Upward_step_power_spectrum_k_c_to_k_s(Upward_step_spectra_k_s, Upward_step_spectra_k_c,prec); //k_s 由 1E-2*k_s 和 k_c 间的 k^4 幂率来定，远小于为10^-2
    
    arb_log(Upward_step_spectra_ln_k_s,Upward_step_spectra_k_s,prec);
    arb_log(Upward_step_spectra_ln_k_c,Upward_step_spectra_k_c,prec);
    arb_inv(Upward_step_spectra_tau_c,Upward_step_spectra_k_c,prec); //τ_c=1/k_c
    
}
 
