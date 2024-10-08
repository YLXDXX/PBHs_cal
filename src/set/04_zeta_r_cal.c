#include "header/04_zeta_r_cal.h"

void Set_zeta_r_cal(slong prec)
{
    //求解 ζ(r) 
    //Ln_K_star=30.37829203018403957048
    //K_star=1.56E13
    //设置求ζ_G(r)辅助变量的积分区域
    switch(Power_spectrum_type)
    {
        case lognormal_type :
            //以ln(k)作为积分变量，若出现负值，正常
            arb_mul_ui(Log_normal_mul_sigma,Power_sigma, 13, prec); //这里的范围需要大点
            arb_sub(Int_sigma_n_min,Ln_K_star,Log_normal_mul_sigma,prec); //动态积分范围
            arb_add(Int_sigma_n_max,Ln_K_star,Log_normal_mul_sigma,prec);
            
            arb_set_str(Int_sigma_n_precision,"1E-12",prec);
            
            break;
            
        case box_type :
            //这里是以ln(k)作为积分变量的
            arb_t t_aa,t_bb,t_cc;
            arb_init(t_aa);
            arb_init(t_bb);
            arb_init(t_cc);
            
            arb_div_ui(t_cc,Box_width,2,prec);
            
            arb_sub(t_aa,Box_K_middle,t_cc,prec); //aa=K_middle - width/2
            arb_log(Int_sigma_n_min,t_aa,prec); //ln(k)作变量，取对数
            
            arb_add(t_bb,Box_K_middle,t_cc,prec); //bb=K_middle + width/2
            arb_log(Int_sigma_n_max,t_bb,prec); //ln(k)作变量，取对数
            
            
            arb_set_str(Int_sigma_n_precision,"1E-8",prec);
            
            break;
            
        case power_law_type :
            //这里是以ln(k)作为积分变量的
            arb_set_str(Int_sigma_n_min,"0",prec);
            arb_set_str(Int_sigma_n_max,"9",prec); //ln(10^13)≈30
            
            arb_set_str(Int_sigma_n_precision,"1E-8",prec);
            
            break;
            
        case broken_power_law_type :
            //以ln(k)作为积分变量，若出现负值，正常
            //在左侧大约以 α 斜率上升，右侧大约经 β 斜率下降，跟据 α/β 值动态取值
            //e^(-n)*k^α=K_star^α => lnk = ln(k_star) -n/α
            arb_ui_div(Int_sigma_n_min,28,BPL_alpha,prec); //左边，n=28
            arb_sub(Int_sigma_n_min,Ln_K_star,Int_sigma_n_min,prec);
            
            arb_ui_div(Int_sigma_n_max,28,BPL_beta,prec); //右边，n=28
            arb_add(Int_sigma_n_max,Ln_K_star,Int_sigma_n_max,prec);
            
            arb_set_str(Int_sigma_n_precision,"1E-12",prec);
            
            break;
        case link_cmb_type :
            //因为前面 Link_CMB_P_m=0, 这里积分区间为[0,k_m]
            //k的最小值可取零，对应ln(k)为-∞，取k=5E-2作截断
            arb_set_str(Int_sigma_n_min,"5E-2",prec);
            arb_log(Int_sigma_n_min,Int_sigma_n_min,prec);//ln(5E-2)
            //arb_neg_inf(Int_sigma_n_min); //直接用负无穷，需要非常高的精度才能直接算
            
            
            //arb_set(Int_sigma_n_max,Link_CMB_K_m);
            arb_set_str(Int_sigma_n_max,"4.6E15",prec); //利用PBHs蒸发对应的质量下限，得到的k_max
            arb_log(Int_sigma_n_max,Int_sigma_n_max,prec);//以ln(k)作变量
            
            
            arb_set_str(Int_sigma_n_precision,"1E-15",prec);
            
            break;
        case numerical_cal_type :
            //以ln(k)作为积分变量，若出现负值，正常
            //这里，积分的上下限由拟合的上下限来决定
            arb_set(Int_sigma_n_min,FITTED_x); //下限
            arb_set(Int_sigma_n_max,FITTED_x+FITTED_num-1); //上限
            arb_set_str(Int_sigma_n_precision,"1E-5",prec); //精度
            
            break;
        case upward_step_spectra_type :
            //并不做PBHs相关的具体计算，仅计算SIGWs，随意设
            arb_set_str(Int_sigma_n_min, "0", prec);
            arb_set_str(Int_sigma_n_max, "1", prec);
            arb_set_str(Int_sigma_n_precision,"1E-12",prec);
            
            break;
        case delta_type :
            //不需要
            break;
        default :
            printf(" main.c Power_spectrum_type 有误\n");
            exit(1);
    }
    
    
    
    //计算ζ(r)相关辅助变量
    //注意，下面的量
    //在PS中，delta情况不需要计算，但在PT中，δ情况下是需要计算的
    //计算(σ_n)^2
    // (σ_n)^2 的量级在 [10*K_star]^n
    
    if(Stdout_verbose==true)
    {
        printf("计算： (σ_0)^2  ---  (σ_4)^2\n");
    }
    
    Help_sigma_n_square(Sigma_0_square,0,prec);
    
    if(Stdout_verbose==true)
    {
        printf("(σ_0)^2 计算完成\n");
    }
    
    Help_sigma_n_square(Sigma_1_square,1,prec);
    
    if(Stdout_verbose==true)
    {
        printf("(σ_1)^2 计算完成\n");
    }
    
    Help_sigma_n_square(Sigma_2_square,2,prec);
    
    if(Stdout_verbose==true)
    {
        printf("(σ_2)^2 计算完成\n");
    }
    
    
    Help_sigma_n_square(Sigma_3_square,3,prec);
    
    if(Stdout_verbose==true)
    {
        printf("(σ_3)^2 计算完成\n");
    }
    
    Help_sigma_n_square(Sigma_4_square,4,prec);
    
    if(Stdout_verbose==true)
    {
        printf("(σ_4)^2 计算完成\n");
    }
    
    
    //计算利用前面计算得到的 σ_n^2 得到 γ_1, γ_3
    //γ_n=σ_n^2/(σ_n-1 * σ_n+1)
    
    arb_t temp_a,temp_b;
    arb_init(temp_a);
    arb_init(temp_b);
    
    // γ_n 的大小在 0.5 左右，基本不受 K_star 影响
    arb_sqrt(temp_a,Sigma_0_square,2*prec); //γ_1
    arb_sqrt(temp_b,Sigma_2_square,2*prec);
    arb_mul(Gamma_1,temp_a,temp_b,prec);
    arb_div(Gamma_1,Sigma_1_square,Gamma_1,prec);
    
    arb_sqrt(temp_a,Sigma_2_square,2*prec); //γ_3
    arb_sqrt(temp_b,Sigma_4_square,2*prec);
    arb_mul(Gamma_3,temp_a,temp_b,prec);
    arb_div(Gamma_3,Sigma_3_square,Gamma_3,prec);
    
    if(Stdout_verbose==true)
    {
        printf("γ_1, γ_3     计算完成\n");
    }
    
    //计算利用前面计算得到的 σ_n^2 得到 R_1, R_3
    //R_n=sqrt(3)*σ_n/σ_n+1
    
    //R_n 的量级在 1/K_star 左右，随 K_star 增大而减小
    arb_sqrt_ui(R_1,3,prec); //R_1
    arb_sqrt(temp_a, Sigma_1_square, prec);
    arb_sqrt(temp_b, Sigma_2_square, prec);
    arb_mul(R_1,R_1,temp_a,prec);
    arb_div(R_1,R_1,temp_b,prec);
    
    arb_sqrt_ui(R_3,3,prec); //R_3
    arb_sqrt(temp_a, Sigma_3_square, prec);
    arb_sqrt(temp_b, Sigma_4_square, prec);
    arb_mul(R_3,R_3,temp_a,prec);
    arb_div(R_3,R_3,temp_b,prec);
    
    if(Stdout_verbose==true)
    {
        printf("R_1, R_3     计算完成\n");
    }
    
    
    arb_clear(temp_a);
    arb_clear(temp_b);
    
    
}
 
