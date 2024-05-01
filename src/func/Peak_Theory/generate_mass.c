#include "generate_mass.h"
#include <stdlib.h>

//在视界重进入时，对应的视界质量 M_H(k)
int Horizon_reentry_M_H(arb_t res,const arb_t k,slong prec)
{
    //函数中所用变量
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    //单位克
    //M_H(k)=10^20 * (g_star/106.75)^(-1/6) * (k / 1.56E13)^(-2)
    //前面系数
    arb_set_str(s,"1E20",prec);
    arb_set(res,s);
    
    //中间
    arb_set_str(s,"106.75",prec);
    arb_div(t,effective_g_star,s,prec);
    arb_one(s);
    arb_div_si(s,s,6,prec);
    arb_neg(s,s);
    arb_pow(t,t,s,prec);
    arb_mul(res,res,t,prec);
    
    //最后
    arb_set_str(s,"1.56E13",prec);
    arb_div(s,k,s,prec);
    arb_sqr(s,s,prec);
    arb_inv(s,s,prec);
    arb_mul(res,res,s,prec);
    
    arb_clear(s);
    arb_clear(t);
    
    return 0;
}


//求质量M，两种表示，一种是相对质量 M/M_H，一种是实际质量 M(μ)
//在视界重进入时，形成的黑洞的质量为M(µ)
int Horizon_reentry_M_mu(arb_t res,const arb_t mu,slong prec)
{
    //函数中所用变量
    arb_t s,t,w;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //
    //修改 Mu_2，非高斯时，r_m 与 Mu_2 有关
    //
    arb_set(Mu_2,mu);
    
    
    //在 M(µ) 近似的可写成如下关系式
    // M(μ) = x_m^2 * exp(2*zeta(mu)) * K * (mu-mu_th)^γ * M_k(K_star)
    
    //最前面系数
    Find_r_max(s, prec); //需重新计算最大值
    arb_set(R_MAX,s); // 赋值，供后面计算 dln(M)/dμ 用
    
    arb_mul(t,s,K_star,prec);
    arb_sqr(w,t,prec);
    
    //指数
    zeta_profile_n(t,s,0,prec);
    arb_mul_si(t,t,2,prec);
    arb_exp(t,t,prec);
    arb_mul(w,w,t,prec);
    
    //中间
    arb_mul(w,w,Mass_K,prec);
    arb_sub(t,mu,Mu_2_th,prec);
    arb_abs(t,t);
    arb_pow(t,t,Mass_gamma,prec);
    arb_mul(w,w,t,prec);
    
    //最后
    Horizon_reentry_M_H(t,K_star,prec);
    arb_mul(res,w,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    
    return 0;
}



//求质量的M，两种表示，一种是相对质量 M/M_H，一种是实际质量 M(μ)
//在视界进入时，视界质量M_H，形成的黑洞的质量为M，两者之比 M/M_H
int Horizon_reentry_M_ratio(arb_t res, const arb_t mu, slong prec)
{
    //函数中所用变量
    arb_t s,r_m,x;
    
    arb_init(s);
    arb_init(r_m);
    arb_init(x);
    
    //修改 Mu_2，非高斯时，最大值与Mu_2有关
    arb_set(Mu_2,mu);
    
    
    //M=x^2 * exp(2*zeta(mu)) * K * (mu-mu_th)^{γ} * M_H
    //前面系数部分
    Find_r_max(r_m, prec); //需重新寻找最大值
    arb_set(R_MAX,r_m); // 赋值，供后面计算 dln(M)/dμ 用
    
    arb_mul(x,r_m,K_star,prec);
    arb_sqr(x,x,prec);
    
    zeta_profile_n(s,r_m,0,prec);
    arb_mul_si(s,s,2,prec);
    arb_exp(s,s,prec);
    
    arb_mul(x,x,s,prec);
    arb_mul(x,x,Mass_K,prec);
    
    
    //幂律部分
    arb_sub(s,mu,Mu_2_th,prec);
    
    arb_abs(s,s); //取绝对值，防止后面开幂次里面为负
    arb_pow(s,s,Mass_gamma,prec);
    arb_mul(res,x,s,prec);
    
    
    arb_clear(s);
    arb_clear(r_m);
    arb_clear(x);
    
    return 0;
}



//在视界重进入时，形成的黑洞的质量为M(µ)的反函数µ(M)
int Horizon_reentry_mu_M(arb_t res,const arb_t M,slong prec)
{
    //修改 Mu_2，非高斯时，最大值与Mu_2有关
    //arb_set(Mu_2,mu);
    
    //对于非高斯的情况，需要特特处理
    //  1，r_max与µ值有关
    //  2，µ(M) 可能不是单值函数，一个M，对应多个µ值
    //对于非高斯的情况应可用M(µ)=0的方式，通过求根求µ(M)
    
    //函数中所用变量
    arb_t s,t,find_M;
    arb_init(s);
    arb_init(t);
    arb_init(find_M);
    
    
    
    //这里采用求解方程零点的算法，来进行求解
    //速度很慢
    
    arb_set(find_M,M);//这里，通参参数 find_M 将 M 的值传入定义函数内部
    
    int zeo_mu_M_func(arb_t res, const arb_t mu, void * t_M, const slong order, slong prec)
    {
        arb_t s;
        arb_init(s);
        
        printf("\n当前 μ_2 值： ");arb_printn(mu, 15,0);
        printf("\t M 值： ");arb_printn(t_M, 15,0);printf("\n计算对应的 M：");
        
        //通过 μ 值计算相应的 M
        if( Relative_Mass )
        {
            Horizon_reentry_M_ratio(s,mu,prec); //使用相对质量的版本
        }else
        {
            Horizon_reentry_M_mu(s,mu,prec); 
        }
        
        arb_printn(s, 15,0);printf("\n两都相差为：  ");
        
        arb_sub(res,s,t_M,prec);
        
        arb_printn(res, 15,0);printf("\n\n");
        
        arb_clear(s);
        return 0;
    }
    
    //arb_set_str(find_min,"0.61532877",prec);
    //由于mu的值需要大于 Mu_2_th ,故利用 Mu_2_th 来设置 mu 的最小值 find_min
    //当M/M_k_*=[0.01, 10]这个区间外时，想要求出对应的mu值
    // - 对于小的 M/M_k 需要减小 find_step
    // - 对于大的 M/M_k 需要增加 find_max
    
    // 这时 Mu_2_th 会随 k_3 变化，每次重新设定一下
    arb_set(Root_M_to_mu_min,Mu_2_th); 
    arb_set_str(t,"1E-25",prec);
    arb_add(Root_M_to_mu_min, Root_M_to_mu_min, t, prec); //精确相等，会导致开根错误
    
    printf("\n\nM -> μ_2 各参数\nM:       ");
    arb_printn(find_M, 15,0);printf("\nMu_2_th: ");
    arb_printn(Mu_2_th, 15,0);printf("\n\n");
    
    Find_interval_root(res, zeo_mu_M_func, find_M, 0,
                       Root_M_to_mu_min, Root_M_to_mu_max, Root_M_to_mu_precision,
                       Root_M_to_mu_num, Root_Normal, prec);
    
    printf("找到 M 所需的 μ_2 ： ");arb_printn(res, 15,0);printf("\n\n");
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(find_M);
    
    return 0;
}



//形成黑洞质M对mu的导数 dln(M)/dµ
int Horizon_reentry_D_ln_M_to_mu(arb_t res, const arb_t mu, slong prec)
{
    //函数中所用变量
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    //跟功率谱有关，对于delta型式的功率谱，可解析求解
    //功率谱类型判断
    switch(Power_spectrum_type) 
    {
        case lognormal_type :
            
            if( SIMPLIFY )
            {
                //简化的情况，其中：
                // r_m 与 k_3 无关，只与 Mu_2 有关
                // 并且，利用 Mu_2 解出 r_m 后，且认为 r_m 是个常数，不参与求导
                // Mu_2_th 与 k_3 无关，（与 Mu_2 无关）
                // 而 r_m 的值，在上一步 Horizon_reentry_mu_M 中
                // 调用 Horizon_reentry_M_mu 或 Horizon_reentry_M_ratio 已解出
                // K(K_3_square) 被认为等于一 K(K_3_square)=1
                // 且在简化的情况下，前面作了如下假设 arb_set(K_3_square,Gamma_3);
                // 此时有 ζ(r)=μ_2 * ψ_1(r) Help_psi_1_n(phi_1,r,0,prec);
                // 并且 μ_2th 与 K_3_square 已无关系
                
                // dln(M)/dµ 与 k_3 无关
                
                // 2*ψ_1(r) + γ/(μ_2-μ_2th)
                
                Help_psi_1_n(s,R_MAX,0,prec); //前面已解出 R_MAX
                arb_mul_si(s,s,2,prec);
                
                arb_sub(t,mu,Mu_2_th,prec);
                arb_div(t,Mass_gamma,t,prec);
                
                arb_add(s,s,t,prec);
                arb_abs(res,s);
                
            } else
            {
                printf("Peak_Theory -> generate_mass -> Horizon_reentry_D_ln_M_to_mu -> lognormal_type 中非近似形式，还未完成\n");
                exit(1);
            }
            
            break;
        case delta_type :
            // 2*sinc(x) + γ/(μ_2-μ_2th) //x=r_m*k_star
            //近似写为 0.282+0.36/(mu-mu_th)
            
            arb_mul(s,R_MAX,K_star,prec);//前面
            arb_sinc(s,s,prec);
            arb_mul_si(s,s,2,prec);
            
            arb_sub(t,mu,Mu_2_th,prec);//后面
            arb_div(t,Mass_gamma,t,prec);
            
            arb_add(s,s,t,prec);
            arb_abs(res,s);
            
            break;
            
        default :
            printf("Peak_Theory -> generate_mass -> power_spectrum_type 有误\n");
            exit(1);
    }
    
    arb_clear(s);
    arb_clear(t);
    
    return 0;
}

 
