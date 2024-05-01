#include "ps_abundance_all_k.h" 
#include <stdlib.h>


//考虑连续谱中所有k模式的影响
//以δ函数求和的形式求连续功率谱的β


//设定 k_0 模式对应的 Σ
int interior_delta_k_Sigma(const arb_t k_0, slong prec) //传过来的值需取对数
{
    //对于δ谱而言，在程序开始运行到此时，已经进行了一遍相关计算
    //可以共用的有：阈值、X_m=r_m*k_*
    //而方差的计算，只差个倍数关系：Σ_new = A_new/A_old * Σ_old
    //其中，在最初的运算中，A_old 和 Σ_old 已经知道了
    //对于log-normal谱，当k=K_star时，对应的振幅为A/(sqrt(2π)*Δ)
    
    //注意，这里是以ln(k)作为变量，k -> ln(k)
    
    arb_t t,s,A_new;
    
    arb_init(t);
    arb_init(s);
    arb_init(A_new);
    
    //对所求功率谱进行判断
    //为加快速度，这里仅进行一次判断
    if(Continuum_spectrum_judge_help!=true) //bool类型默认为假
    {
        if((Continuum_spectrum_cal_simplify==true && Power_spectrum_type!=delta_type) ||
            (Continuum_spectrum_cal_simplify!=true && Power_spectrum_type==delta_type) ||
            Continuum_spectrum_type==delta_type)
        {
            printf("连续谱计算，条件冲突：\n");
            printf("\t简化为真，Power_spectrum_type 需要是 delta_type\n");
            printf("\t简化为假，Power_spectrum_type 不能为 delta_type\n");
            printf("\t且 Continuum_spectrum_type 不能为 delta_type\n");
            exit(1);
        }
        //printf("连续谱计算，条件良好\n");
    }
    
    POWER spectrum_type_help;//交换spectrum_type与Power_spectrum_type
    spectrum_type_help=Continuum_spectrum_type;
    Continuum_spectrum_type=Power_spectrum_type;
    Power_spectrum_type=spectrum_type_help;
    
    //关于这里 A_old 的求解，对称和非对称的谱，求解方式不同
    //为加快速度，这里仅对 A_old 进行一次计算
    if(Continuum_spectrum_judge_help!=true) //bool类型默认为假
    {
        arb_t t_k;
        arb_init(t_k);
        
        switch(Power_spectrum_type) //作了变量交换，使用Power_spectrum_type判断
        {
            case lognormal_type : //取Σ_old对应的模式为x_m/r_m，其中x_m=r_m×k_*是δ情况的取值
                if(Continuum_spectrum_cal_simplify==true)
                {
                    power_spectrum(Continuum_spectrum_A_old, Ln_K_star, prec); //获限A_old, 以 ln(k) 作参数传入
                }else{
                    //arb_sqrt(t,Pi_2,prec); // A_old=A/(sqrt(2π)*Δ)
                    //arb_mul(t,t,Power_sigma,prec);
                    //arb_div(Continuum_spectrum_A_old,Power_A,t,prec);
                    //直接从功率谱中得到，对功率谱的修改可以直接起作用，保持一致性
                    //采取δ形式的k*r_m=const.
                    arb_log(t_k,Continuum_spectrum_k_ch,prec); //取对数
                    power_spectrum(Continuum_spectrum_A_old, t_k, prec); //获限A_old, 以 ln(k) 作参数传入
                }
                break;
            case power_law_type : //取Σ_old对应的模式为x_m/r_m，其中x_m=r_m×k_*是δ情况的取值
                
                arb_log(t_k,Continuum_spectrum_k_ch,prec); //取对数
                power_spectrum(Continuum_spectrum_A_old, t_k, prec); //获限A_old, 以 ln(k) 作参数传入
                
                break;
            case box_type : //取Σ_old对应的模式为x_m/r_m，其中x_m=r_m×k_*是δ情况的取值
                
                //arb_set(Continuum_spectrum_A_old,Power_A);// A_old=Power_A
                //直接从功率谱中得到，对功率谱的修改可以直接起作用，保持一致性
                power_spectrum(Continuum_spectrum_A_old, Box_K_middle, prec); //获限A_old, 以 ln(k) 作参数传入
                
                break;
            case broken_power_law_type : //取Σ_old对应的模式为x_m/r_m，其中x_m=r_m×k_*是δ情况的取值
                
                arb_log(t_k,Continuum_spectrum_k_ch,prec); //取对数
                power_spectrum(Continuum_spectrum_A_old, t_k, prec); //获限A_old, 以 ln(k) 作参数传入
                
                break;
            case link_cmb_type : //取Σ_old对应的模式为x_m/r_m，其中x_m=r_m×k_*是δ情况的取值
                
                arb_log(t_k,Continuum_spectrum_k_ch,prec); //取对数
                power_spectrum(Continuum_spectrum_A_old, t_k, prec); //获限A_old, 以 ln(k) 作参数传入
                
                break;
            default :
                printf("Continuum_spectrum_type 有误，无法正确计算特征模式\n");
                exit(1);
        }
        arb_clear(t_k);
        Continuum_spectrum_judge_help=true;
    }
    
    power_spectrum(A_new, k_0, prec); //获限A_new, 以 ln(k) 作参数传入
    
    spectrum_type_help=Continuum_spectrum_type;//再次交换，复位
    Continuum_spectrum_type=Power_spectrum_type;
    Power_spectrum_type=spectrum_type_help;
    
    //对于有些功率谱，在某范围之外，振幅为零
    if(arb_is_zero(A_new))
    {
        arb_clear(t);
        arb_clear(s);
        arb_clear(A_new);
        
        return 1;
    }
    
    arb_div(t,A_new,Continuum_spectrum_A_old,prec); //Σ_new = A_new/A_old * Σ_old
    arb_mul(PS_Sigma_XX,PS_Sigma_XX_save,t,prec);
    arb_mul(PS_Sigma_XY,PS_Sigma_XY_save,t,prec);
    arb_set(PS_Sigma_YX,PS_Sigma_XY);
    arb_mul(PS_Sigma_YY,PS_Sigma_YY_save,t,prec);
    
    
    arb_clear(t);
    arb_clear(s);
    arb_clear(A_new);
    
    return 0;
}

//求单个 k_0 模式下总的 β，β=∫β(M/M_H) d M/M_H
int PS_abundance_beta_delta_k(arb_t res, const arb_t k_0, slong prec) //传过来的值需取对数
{
    
    arb_t t;
    arb_init(t);
    
    
    int i;
    i=interior_delta_k_Sigma(k_0,prec); //设定连续谱中k_0模式协方差，以 ln(k_0) 传入
    
    //判断返回值
    if(i==0) //成功设定k_0模式协方差
    {
        PS_abundance_beta_all(t,prec);
        arb_set(res,t);
    }else{
        arb_zero(res); //设结果为零
    }
    
    
    arb_clear(t);
    
    return 0;
}


int interior_PS_abundance_beta_delta_k_all(arb_t res, const arb_t k_0, void* p, const slong order, slong prec)
{
    arb_t t,s;
    
    arb_init(t);
    arb_init(s);
    
    //注意，这里是以ln(k)作为变量，k -> ln(k)
    //以 k_0 作为积分变量
    
    //∫β(k) dk = ∫β(ln(k))*k*dk/k = ∫β(ln(k))*e^ln(k)*dln(k)
    
    //arb_log(t,k_0,prec);
    PS_abundance_beta_delta_k(t,k_0,prec); //以对数值传入 ln(k_0)
    arb_exp(s,k_0,prec);
    
    arb_mul(res,t,s,prec);
    
    arb_clear(t);
    arb_clear(s);
    
    return 0;
}



int PS_abundance_beta_delta_k_all(arb_t res, slong prec)
{
    arb_t t,s;
    
    arb_init(t);
    arb_init(s);
    
    
    //使用新的gauss_kronrod积分算法
    Integration_arb(res, interior_PS_abundance_beta_delta_k_all, NULL, 0, 
                    PS_abundance_beta_delta_k_all_min, PS_abundance_beta_delta_k_all_max,
                    PS_abundance_beta_delta_k_all_precision,
                    PS_abundance_beta_delta_k_all_Int_iterate_min,PS_abundance_beta_delta_k_all_Int_iterate_max, prec);
    
    
    arb_clear(t);
    arb_clear(s);
    
    return 0;
}


//求某个质量 M 在 PBHs 生成时的占比，考虑所有k模式
//实际上是用质量 M 对应的特征模式 k_M 来进行计算

int interior_PS_abundance_beta_delta_k_M(arb_t res, const arb_t k, void* parameter, const slong order, slong prec)
{
    arb_t t,s;
    
    arb_init(t);
    arb_init(s);
    
    //需传入结构体 Delta_Beta_M_All_k 来获取参数 k_M
    struct Delta_Beta_M_All_k *Beta_M_parameter; //这里不需要手动分配，只需将其指向传入的指针即可
    Beta_M_parameter=parameter;
    
    
    //β[(k/k_M)^3]，由于传入的是取对数后的值
    //故有：(k/k_M)^3=exp( 3*[ln(k)-ln(k_M)] )
    
    arb_sub(t,k,Beta_M_parameter->k_M,prec); 
    arb_mul_ui(t,t,3,prec);
    arb_exp(t,t,prec);
    
    int i;
    i=interior_delta_k_Sigma(k,prec); //设定连续谱中k模式协方差，以对数值传入 ln(k)
    
    //判断返回值
    if(i==0) //成功设定k模式协方差
    {
        PS_abundance_beta_m(s, t, prec); //β[(k/k_M)^3]
        arb_set(res,s);
    }else{
        arb_zero(res); //设结果为零
    }
    
    
    arb_clear(t);
    arb_clear(s);
    
    return 0;
}

int PS_abundance_beta_delta_k_M(arb_t res, const arb_t k_M, slong prec) //传过来的值需取对数
{
    arb_t t,s,k_up,k_down;
    
    arb_init(t);
    arb_init(s);
    arb_init(k_up);
    arb_init(k_down);
    
    //求积分 k 的上限 k_up= k_M * (M/M_H)^{1/3}  // 其中M/M_H取最大值
    //但传入的是k_M的对数值，ln(k_M)
    //故有：ln(k_up)=ln(k_M) + 1/3*ln(M/M_H)
    
    arb_one(t);  // 1/3*ln(M/M_H)
    arb_div_ui(t,t,3,prec);
    arb_log(s,PS_M_ratio_max,prec);
    arb_mul(s,s,t,prec);
    
    arb_add(k_up,s,k_M,prec); //ln(k_M) + 1/3*ln(M/M_H)
    
    
    //设定积分下限
    switch(Continuum_spectrum_type) //作了变量交换，使用Power_spectrum_type判断
    {
        case lognormal_type:
            
            //对log-normal谱的上界进行优化
            //当积分上下界与k_*差太多，不对称时，double_exponential积分的计算会出问题，需要采用gauss_kronrod_iterat积分
            //但是，上界与中心值差太多，可以进行截断，其中取9σ已足够精确
            
            arb_mul_ui(t, Power_sigma, 9, prec);
            arb_add(t,t,Ln_K_star,prec);
            if(arb_gt(k_up,t)) //若 k_up > Ln_K_star + 9σ 则进行截断
            {
                arb_set(k_up,t);
            }
            
            //对log-normal谱的下界进行优化
            if(arb_gt(k_M,Ln_K_star))
            {
                arb_mul_ui(t,Power_sigma, 9, prec); //k_down=Ln_K_star-9*Δ
                arb_sub(k_down,Ln_K_star,t,prec); 
            }else
            {
                arb_mul_ui(t,Power_sigma, 9, prec); //k_down=k_M-9*Δ
                arb_sub(k_down,k_M,t,prec);
                
                
                //对log-normal谱的下界进行优化
                arb_sub(t,Ln_K_star,t,prec); 
                if( arb_lt(k_down,t) && arb_gt(k_up,t) )//若 k_down < Ln_K_star - 9σ 并且 k_up > t  则进行截断
                {
                    arb_set(k_down,t);
                }
                
            }
            break;
        default :
            printf("Continuum_spectrum_type 有误，当前仅对lognormal_type进行了适配\n");
            exit(1);
    }
    
    
    //利用结构体 Delta_Beta_M_All_k 向被积函数传参
    //需手动分配内存
    struct Delta_Beta_M_All_k *Beta_M_parameter = (struct Delta_Beta_M_All_k *)calloc(1,sizeof(struct Delta_Beta_M_All_k));
    
    arb_init(Beta_M_parameter->k_M);//使用arb_t变量前初始化
    
    arb_set(Beta_M_parameter->k_M,k_M); //设定求根参数
    
    //使用新的gauss_kronrod积分算法
    Integration_arb(t, interior_PS_abundance_beta_delta_k_M, Beta_M_parameter, 0, 
                    k_down, k_up,
                    PS_abundance_beta_delta_k_M_precision,
                    PS_abundance_beta_delta_k_M_Int_iterate_min,PS_abundance_beta_delta_k_M_Int_iterate_max, prec);
    
    arb_set(res,t);
    
    arb_clear(t);
    arb_clear(s);
    arb_clear(k_up);
    arb_clear(k_down);
    
    return 0;
}
