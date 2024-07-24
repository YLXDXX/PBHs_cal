#include "pt_generate_mass.h"
#include <stdlib.h>


struct M_TO_MU {
  arb_t zeta_k;
  arb_t M;
};


//k模式进入视界，所对应的质量 M_(k)
int Horizon_reentry_k_to_M(arb_t res,const arb_t k,slong prec)
{
    //函数中所用变量
    arb_t s,t,w;
    
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //注意到，单位为 g （克），2109.00791 (3.7)
    //M_H(k)=10^20 * (g_star/106.75)^(-1/6) * (k / 1.56E13)^(-2)
    
    arb_set_str(w,"1E20",prec); //前面系数
    
    Func_k_to_degrees_of_freedom(t, s, k, prec); //求k对应的自由度数,第一个是自由度数，第二个是熵自由度数
    
    arb_set_str(s,"106.75",prec); //中间
    arb_div(t,t,s,prec);
    
    arb_one(s);
    arb_div_si(s,s,6,prec);
    arb_neg(s,s);
    arb_pow(t,t,s,prec);
    arb_mul(w,w,t,prec);
    
    
    arb_set_str(s,"1.56E13",prec); //最后
    arb_div(s,k,s,prec);
    arb_sqr(s,s,prec);
    arb_inv(s,s,prec);
    
    arb_mul(res,w,s,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    
    return 0;
}


//求质量M，两种表示，一种是相对质量 M/M_H，一种是实际质量 M(μ)
//在视界重进入时，形成的黑洞的质量为M(µ)
int Horizon_reentry_mu_to_M(arb_t res, const arb_t mu, const arb_t zeta_k, slong prec)
{
    //函数中所用变量
    arb_t s,t,w,mu_th,r_m;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(mu_th);
    arb_init(r_m);
    
    //修改 PT_mu，非高斯时，r_m 与 PT_mu 有关
    //profile非简化时，r_m 还与 k 有关
    
    arb_set(PT_mu,mu);
    arb_set(PT_k,zeta_k);
    
    if(PT_threshold_simplify==true)
    {
       arb_set(mu_th,PT_mu_th); 
    }else
    {
        Find_PT_Mu_th(mu_th, zeta_k, prec); //计算 μ_th 供后面使用
        arb_set(PT_mu_th,mu_th);
    }
    
    
    //在 M(µ) 近似的可写成如下关系式 2109.00791 (3.8)
    // M(μ) = x_m^2 * exp(2*zeta(mu)) * K * (C-C_th)^γ * M_k(K_star)
    
    //最前面系数
    if(PT_cal_r_m_fix==true)
    {
        //虽然 R_MAX 与 μ和k都有关，不再重新计算 R_MAX
        arb_set(r_m,R_MAX);
    }else
    {
        Find_r_max(r_m, prec); //需重新计算最大值
        arb_set(R_MAX,r_m); // 赋值，供后面计算 dln(M)/dμ 用
    }
    
    
    arb_mul(t,r_m,K_star,prec);
    arb_sqr(w,t,prec);
    
    //指数
    zeta_profile_n(t,r_m,0,prec);
    arb_mul_si(t,t,2,prec);
    arb_exp(t,t,prec);
    arb_mul(w,w,t,prec);
    
    //中间
    arb_mul(w,w,Mass_K,prec);
    
    //中间，幂律部分
    //原论文中使用的是 μ-μ_c，已改成 C-C_th
    C_r_profile_n(s,r_m,0,prec);
    
    arb_set(PT_mu,mu_th); //注意，这里不能直接使用 C_th
                        //因为求解 C_th 时并没有直接使用 C(mu_th) 两都有一个 1E-10 的差
                        //若使用 C_th 对于较小的M不能求出对应的μ
                        //另外，对于 mu_th 也可能变动
    C_r_profile_n(t,r_m,0,prec);
    arb_set(PT_mu,mu); //求解完后，重新设为当前的 μ
    
    arb_sub(s,s,t,prec);
    
    arb_abs(s,s);
    arb_pow(s,s,Mass_gamma,prec);
    arb_mul(w,w,s,prec);
    
    //最后
    //注意到，M_H 与r_m相关，而r_m又与μ_2和k_3相关
    Horizon_reentry_k_to_M(t,K_star,prec);
    arb_mul(res,w,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(mu_th);
    arb_clear(r_m);
    
    return 0;
}



//求质量的M，两种表示，一种是相对质量 M/M_H，一种是实际质量 M(μ)
//在视界进入时，视界质量M_H，形成的黑洞的质量为M，两者之比 M/M_H
int Horizon_reentry_mu_to_M_relative(arb_t res, const arb_t mu, const arb_t zeta_k, slong prec)
{
    //函数中所用变量
    arb_t s,t,r_m,x,mu_th;
    
    arb_init(s);
    arb_init(t);
    arb_init(r_m);
    arb_init(x);
    arb_init(mu_th);
    
    if(PT_threshold_simplify==true)
    {
        arb_set(mu_th,PT_mu_th); 
    }else
    {
        Find_PT_Mu_th(mu_th, zeta_k, prec); //计算 μ_th 供后面使用
        arb_set(PT_mu_th,mu_th);
    }
    
    //修改 PT_mu，非高斯时，r_m 与 PT_mu 有关
    //profile非简化时，r_m 还与 k 有关
    
    arb_set(PT_mu,mu);
    arb_set(PT_k,zeta_k);
    
    
    //M=x^2 * exp(2*zeta(mu)) * K * (C-C_th)^{γ} * M_k(k_star)
    //前面系数部分
    if(PT_cal_r_m_fix==true)
    {
        //虽然 R_MAX 与 μ和k都有关，不再重新计算 R_MAX
        arb_set(r_m,R_MAX);
    }else
    {
        Find_r_max(r_m, prec); //需重新计算最大值
        arb_set(R_MAX,r_m); // 赋值，供后面计算 dln(M)/dμ 用
    }
    
    arb_mul(x,r_m,K_star,prec);
    arb_sqr(x,x,prec);
    
    zeta_profile_n(s,r_m,0,prec);
    arb_mul_si(s,s,2,prec);
    arb_exp(s,s,prec);
    
    arb_mul(x,x,s,prec);
    arb_mul(x,x,Mass_K,prec);
    
    
    //幂律部分
    //原论文中使用的是 μ-μ_c，已改成 C-C_th
    C_r_profile_n(s,r_m,0,prec);
    
    arb_set(PT_mu,mu_th); //注意，这里不能直接使用 C_th
                        //因为求解 C_th 时并没有直接使用 C(mu_th) 两都有一个 1E-10 的差
                        //若使用 C_th 对于较小的M不能求出对应的μ
                        //另外，对于 mu_th 也可能变动
    C_r_profile_n(t,r_m,0,prec);
    arb_set(PT_mu,mu); //求解完后，重新设为当前的 μ
    
    
    arb_sub(s,s,t,prec);
    
    
    arb_abs(s,s); //取绝对值，防止后面开幂次里面为负
    arb_pow(s,s,Mass_gamma,prec);
    arb_mul(res,x,s,prec);
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(r_m);
    arb_clear(x);
    arb_clear(mu_th);
    
    return 0;
}

//通过M反解出μ
static int interior_M_to_mu_func(arb_t res, const arb_t mu, void * param, const slong order, slong prec)
{
    arb_t s;
    arb_init(s);
    
    struct M_TO_MU *parameter;
    
    parameter=param;
    
    if(Stdout_verbose==true)
    {
        printf("\n当前 μ_2 值： ");arb_printn(mu, 15,0);
        printf("\t M 值： ");arb_printn(parameter->M, 15,0);printf("\n计算对应的 M："); 
    }
    
    
    //通过 μ 值计算相应的 M
    if( PT_Mass_Relative )
    {
        Horizon_reentry_mu_to_M_relative(s,mu,parameter->zeta_k,prec); //使用相对质量的版本
    }else
    {
        Horizon_reentry_mu_to_M(s,mu,parameter->zeta_k,prec); 
    }
    
    arb_sub(res,s,parameter->M,prec);
    
    if(Stdout_verbose==true)
    {
        arb_printn(s, 15,0);printf("\n两都相差为：  ");
        arb_printn(res, 15,0);printf("\n\n");
    }
    
    arb_clear(s);
    return 0;
}

//在视界重进入时，形成的黑洞的质量为M(µ)的反函数µ(M)
int Horizon_reentry_M_to_mu(arb_t res, const arb_t M, const arb_t zeta_k, slong prec)
{
    //对于非高斯的情况，需要特特处理
    //  1，r_max与µ值有关
    //  2，µ(M) 可能不是单值函数，一个M，对应多个µ值
    //对于非高斯的情况应可用M(µ)=0的方式，通过求根求µ(M)
    
    //这里，没有近似的解析表达式可用
    //用求根方法求解，耗时很长，注意与之相关积分的参数设定
    
    //函数中所用变量
    arb_t s,t,a,b,ratio,r_m;
    arb_init(s);
    arb_init(t);
    arb_init(a);
    arb_init(b);
    arb_init(r_m);
    arb_init(ratio);
    
    
    //通参结构体指针 将 M,k 的值传入定义函数内部
    //这里，对于结构体 Find_root_delta_C_l_Y 需手动分配内存
    struct M_TO_MU *parameter = (struct M_TO_MU *)calloc(1,sizeof(struct M_TO_MU));
    
    arb_init(parameter->zeta_k);//使用arb_t变量前初始化
    arb_init(parameter->M);
    
    arb_set(parameter->zeta_k,zeta_k);
    arb_set(parameter->M,M);
    
    //这里采用求解方程零点的算法，来进行求解
    //arb_set_str(find_min,"0.61532877",prec);
    //由于mu的值需要大于 PT_mu_th ,故利用 PT_mu_th 来设置 mu 的最小值 find_min
    //当M/M_k_*=[0.01, 10]这个区间外时，想要求出对应的mu值
    // - 对于小的 M/M_k 需要减小 find_step
    // - 对于大的 M/M_k 需要增加 find_max
    
    // 这时 PT_mu_th 会随 k 变化，每次重新设定一下
    arb_set(a,PT_mu_th); 
    arb_set_str(t,"1E-40",prec);
    arb_add(a, a, t, prec); //精确相等，会导致开根错误
    
    arb_set(b,PT_mu_max);
    arb_set_str(t,"0.5",prec);//这里多一点，测试不同M时用
    arb_add(b,b,t,prec);
    
    if(Stdout_verbose==true)
    {
        printf("\n\nM -> μ_2 各参数\nM:       ");
        arb_printn(parameter->M, 15,0);printf("\nPT_mu_th: ");
        arb_printn(PT_mu_th, 15,0);printf("\n\n");
    }
    
    
    Find_interval_root(s, interior_M_to_mu_func, parameter, 0,
                       a, b, Root_M_to_mu_precision,
                       Root_M_to_mu_num, Root_Normal, prec);
    arb_set(res,s);
    
    if(Stdout_verbose==true)
    {
        printf("找到 M 所需的 μ_2 ： ");arb_printn(res, 15,0);printf("\n\n");
    } 
    
    //arb_printn(res, 50,0);printf("\n");
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(a);
    arb_clear(b);
    arb_clear(r_m);
    arb_clear(ratio);
    
    return 0;
}



//形成黑洞质M对mu的导数 dln(M)/dµ
//注意，这里的求导，相对质量和非相对质量，所对应的结果都是一样的
int Horizon_reentry_derivative_ln_M_mu(arb_t res, const arb_t mu, const arb_t zeta_k, slong prec)
{
    //函数中所用变量
    arb_t s,t,q,w,C,C_th,C_l,r_m,mu_th,zeta_G_r,zeta_G_r_prime,cal_c;
    arb_init(s);
    arb_init(t);
    arb_init(q);
    arb_init(w);
    arb_init(C);
    arb_init(C_th);
    arb_init(C_l);
    arb_init(r_m);
    arb_init(mu_th);
    arb_init(zeta_G_r);
    arb_init(zeta_G_r_prime);
    arb_init(cal_c);
    
    int label;
    label=0;
    
    arb_set(PT_mu,mu); //导数的值与 μ 相关
    arb_set(PT_k,zeta_k); //导数的值可能跟 k 有关
    
    
    arb_set(r_m,R_MAX); //前面调用的函数已解出 R_MAX，这里可直接用
    arb_set(mu_th,PT_mu_th); //前面调用的函数已解出 μ_th，这里可直接用
    
    //跟功率谱有关，对于delta型式的功率谱，可解析求解
    //功率谱类型判断
    switch(Power_spectrum_type) 
    {
        //注意到，我们这里导数的计算公式 2109.00792 (3.9) 是将 M_H 看作一个常数的结果
        // M_H --> M(k_*) 这在δ谱下，是成立的，但在连续谱下即使采用profile简化
        //  M_H 也会跟 μ_2 有关
        //注意到，这里进入视界用的是 k_*，通常对应于log-normal和BPL谱的中心尺度
        //对于其它的谱，需要选择某个合适的值，k_*选择的影响，在 x_m=r_m k_* 中体现出来
        
        case lognormal_type :
        case broken_power_law_type :
        case power_law_type :
        case box_type :
        case link_cmb_type :
            
            label=11;
            
            break;
            
        case delta_type :
            // dζ/dμ 高斯情况，与k无关，与μ无关
            // dζ/dμ 非高斯下，与k无关，与μ有关
            
            label=22;
            
            break;
            
        default :
            printf("Peak_Theory -> generate_mass -> Horizon_reentry_derivative_ln_M_mu -> power_spectrum_type 错误 \n");
            exit(1);
    }
    
    
    //ζ_G=μ*ζ_G_r -- ζ=F(ζ_G)
    if(label==11) //连续谱
    {
        if( PT_profile_simplify )
        {
            //简化的情况，其中：
            // r_m 与 k 无关，只与 PT_mu 有关
            // 并且，利用 PT_mu 解出 r_m 后，且认为 r_m 是个常数，不参与求导
            // PT_mu_th 与 k 无关，（也与 PT_mu 无关）
            // 而 r_m 的值，在上一步 Horizon_reentry_M_to_mu 中
            // 调用 Horizon_reentry_mu_to_M 或 Horizon_reentry_mu_to_M_relative 已解出
            // K(K_square) 被认为等于一 K(K_square)=1
            // 且在简化的情况下，前面作了如下假设 arb_set(K_square,Gamma_3);
            // 此时有 ζ(r)=μ_2 * ψ_1(r) --> dζ/dμ=ψ_1(r)
            // 并且 μ_th 与 K_square 已无关系
            // dln(M)/dµ 与 k 无关
            
            // 2*ψ_1(r) + γ/(μ_2-μ_2th)
            
            Help_psi_n(zeta_G_r,r_m,0,prec);
            Help_psi_n(zeta_G_r_prime,r_m,1,prec);
            
        } else
        {
            // dζ/dμ 高斯情况，与k有关，与μ无关
            // dζ/dμ 非高斯下，与k有关，有μ有关
            zeta_Gauss_profile_n_div_mu(zeta_G_r,r_m,0,prec);
            zeta_Gauss_profile_n_div_mu(zeta_G_r_prime,r_m,1,prec);
        }
        
    }else //δ谱: ζ_G_r=sinc(x) 而 x=r_m*r_star
    {
        Help_psi_n(zeta_G_r,r_m,0,prec);
        Help_psi_n(zeta_G_r_prime,r_m,1,prec);
    }
    
    
    //高斯和非高斯，公共部分 γ/(C-C_th)*(1-3/4*C_l)*(-4/3*r)
    
    C_r_profile_n(C,r_m,0,prec); // γ/(C-C_th)
    Trans_C_to_C_l(C_l,C,prec); // C --> C_l
    
    arb_set(PT_mu,mu_th); //注意，这里不能直接使用 C_th
    C_r_profile_n(C_th,r_m,0,prec);
    arb_set(PT_mu,mu); //求解完后，重新设为当前的 μ
    
    
    arb_sub(s,C,C_th,prec);
    arb_div(s,Mass_gamma,s,prec);
    
    arb_mul_si(t,C_l,-3,prec); //(1-3/4*C_l)
    arb_div_ui(t,t,4,prec);
    arb_add_ui(t,t,1,prec);
    arb_mul(s,s,t,prec);
    
    arb_mul_si(t,r_m,-4,prec);//(-4/3*r)
    arb_div_ui(t,t,3,prec);
    arb_mul(cal_c,s,t,prec);
    
    switch(Zeta_type) 
    {
        case gaussian_type :
            //dζ/dμ=ζ_G_r
            // 2*dζ/dμ + γ/(C-C_th)*dC/dμ_2
            //近似写为 0.282+0.36/(mu-mu_th)
            
            //前面 2*dζ/dμ 部分
            arb_mul_si(w,zeta_G_r,2,prec);
            
            
            //后面部分：γ/(C-C_th)*(1-3/4*C_l)*(-4/3*r)*∂/∂μ_2*∂/∂r*F(μ*ζ_G_r)
            //∂/∂μ_2*∂/∂r*F(μ*ζ_G_r)=ζ_G_r'
            arb_mul(s,cal_c,zeta_G_r_prime,prec);
            
            arb_add(w,w,s,prec);
            arb_abs(res,w);//行列式绝对值
            
            break;
        case exponential_tail_type :
            // dζ/dμ=ζ_G_r/(1-β*ζ_G_r*μ)
            
            arb_mul(s,Exponential_tail_beta,zeta_G_r,prec); //(1-β*ζ_G_r*μ)
            arb_mul(q,s,mu,prec);
            arb_neg(s,q);
            arb_add_ui(w,s,1,prec);
            
            arb_div(s,zeta_G_r,w,prec); 
            arb_mul_si(s,s,2,prec); // 2*dζ/dμ 部分完
            
            //∂/∂μ_2*∂/∂r*F(μ*ζ_G_r)
            //左边
            arb_div(t,zeta_G_r_prime,w,prec);
            
            //右边
            arb_sqr(w,w,prec); //w用掉
            arb_mul(q,q,zeta_G_r_prime,prec); //q用掉
            arb_div(q,q,w,prec);
            arb_add(t,t,q,prec); //∂/∂μ_2*∂/∂r*F(μ*ζ_G_r) 完
            
            arb_mul(t,cal_c,t,prec);
            
            arb_add(s,s,t,prec);
            arb_abs(res,s);//行列式绝对值
            
            break;
        case up_step_type:
            // dζ/dμ=ζ_G_r/sqrt(1-h*ζ_G_r*μ)
            arb_t h;
            arb_init(h);
            
            arb_abs(h,Up_step_h);
            
            arb_mul(s,h,zeta_G_r,prec); //sqrt(1-h*ζ_G_r*μ)
            arb_mul(q,s,mu,prec);
            arb_neg(s,q);
            arb_add_ui(w,s,1,prec);
            arb_sqrt(t,w,prec);
            
            arb_div(s,zeta_G_r,t,prec);
            arb_mul_si(s,s,2,prec); // 2*dζ/dμ 部分完
            
            //∂/∂μ_2*∂/∂r*F(μ*ζ_G_r)
            //右边
            arb_mul(w,t,w,prec); //w用掉
            arb_mul_ui(w,w,2,prec);
            arb_mul(q,q,zeta_G_r_prime,prec); //q用掉
            arb_div(q,q,w,prec);
            
            //左边
            arb_div(t,zeta_G_r_prime,t,prec);
            
            arb_add(t,t,q,prec);
            arb_mul(t,cal_c,t,prec); //∂/∂μ_2*∂/∂r*F(μ*ζ_G_r) 完
            
            
            arb_add(s,s,t,prec);
            arb_abs(res,s);//行列式绝对值
            
            arb_clear(h);
            break;
        case power_expansion_type:
            // dζ/dμ=6*E*ζ_G_r^6*μ^5+5*D*ζ_G_r^5*μ^4+4*C*ζ_G_r^4*μ^3+3*B*ζ_G_r^3*μ^2+2*A*ζ_G_r^2*μ+ζ_G_r
            arb_t w,A,B,C,D,E;
            
            arb_init(w);
            arb_init(A);
            arb_init(B);
            arb_init(C);
            arb_init(D);
            arb_init(E);
            
            arb_one(A);
            arb_mul_ui(A,A,3,prec);
            arb_div_ui(A,A,5,prec);
            arb_mul(A,A,Power_expansion_f,prec);
            
            arb_one(B);
            arb_mul_ui(B,B,9,prec);
            arb_div_ui(B,B,25,prec);
            arb_mul(B,B,Power_expansion_g,prec);
            
            arb_set(C,Power_expansion_four);
            arb_set(D,Power_expansion_five);
            arb_set(E,Power_expansion_six);
            
            arb_mul_ui(s,E,6,prec);//6*E*ζ_G_r^6*μ^5
            arb_pow_ui(t,zeta_G_r,6,prec);
            arb_mul(s,s,t,prec);
            arb_pow_ui(t,mu,5,prec);
            arb_mul(w,s,t,prec);
            
            arb_mul_ui(s,D,5,prec);//5*D*ζ_G_r^5*μ^4
            arb_pow_ui(t,zeta_G_r,5,prec);
            arb_mul(s,s,t,prec);
            arb_pow_ui(t,mu,4,prec);
            arb_mul(s,s,t,prec);
            arb_add(w,w,s,prec);
            
            arb_mul_ui(s,C,4,prec);//4*C*ζ_G_r^4*μ^3
            arb_pow_ui(t,zeta_G_r,4,prec);
            arb_mul(s,s,t,prec);
            arb_pow_ui(t,mu,3,prec);
            arb_mul(s,s,t,prec);
            arb_add(w,w,s,prec);
            
            arb_mul_ui(s,B,3,prec);//3*B*ζ_G_r^3*μ^2
            arb_pow_ui(t,zeta_G_r,3,prec);
            arb_mul(s,s,t,prec);
            arb_pow_ui(t,mu,2,prec);
            arb_mul(s,s,t,prec);
            arb_add(w,w,s,prec);
            
            arb_mul_ui(s,A,2,prec);//2*A*ζ_G_r^2*μ
            arb_pow_ui(t,zeta_G_r,2,prec);
            arb_mul(s,s,t,prec);
            arb_mul(s,s,mu,prec);
            arb_add(w,w,s,prec);
            
            //ζ_G_r
            arb_add(t,w,zeta_G_r,prec);
            arb_mul_si(q,t,2,prec); // 2*dζ/dμ 部分完
            
            //∂/∂μ_2*∂/∂r*F(μ*ζ_G_r)
            
            arb_mul_ui(s,E,36,prec);//36*E*ζ_G_r^5*μ^5
            arb_pow_ui(t,zeta_G_r,5,prec);
            arb_mul(s,s,t,prec);
            arb_pow_ui(t,mu,5,prec);
            arb_mul(w,s,t,prec);
            
            arb_mul_ui(s,D,25,prec);//25*D*ζ_G_r^4*μ^4
            arb_pow_ui(t,zeta_G_r,4,prec);
            arb_mul(s,s,t,prec);
            arb_pow_ui(t,mu,4,prec);
            arb_mul(s,s,t,prec);
            arb_add(w,w,s,prec);
            
            arb_mul_ui(s,C,16,prec);//16*C*ζ_G_r^3*μ^3
            arb_pow_ui(t,zeta_G_r,3,prec);
            arb_mul(s,s,t,prec);
            arb_pow_ui(t,mu,3,prec);
            arb_mul(s,s,t,prec);
            arb_add(w,w,s,prec);
            
            arb_mul_ui(s,B,9,prec);//9*B*ζ_G_r^2*μ^2
            arb_pow_ui(t,zeta_G_r,2,prec);
            arb_mul(s,s,t,prec);
            arb_pow_ui(t,mu,2,prec);
            arb_mul(s,s,t,prec);
            arb_add(w,w,s,prec);
            
            
            arb_mul_ui(s,A,4,prec);//4*A*ζ_G_r*μ
            arb_mul(s,s,zeta_G_r,prec);
            arb_mul(s,s,mu,prec);
            arb_add(w,w,s,prec);
            
            arb_add_ui(w,w,1,prec);
            
            arb_mul(w,w,zeta_G_r_prime,prec); //最后乘上 ζ_G_r_prime
            
            arb_mul(w,cal_c,w,prec); //∂/∂μ_2*∂/∂r*F(μ*ζ_G_r) 完
            
            arb_add(s,q,w,prec);
            arb_abs(res,s);//行列式绝对值
            
            arb_clear(w);
            arb_clear(A);
            arb_clear(B);
            arb_clear(C);
            arb_clear(D);
            arb_clear(E);
            break;
            default:
                printf("Peak_Theory -> generate_mass -> Horizon_reentry_derivative_ln_M_mu 中 zeta_type 输入有误\n");
                exit(1);
    }
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(q);
    arb_clear(w);
    arb_clear(C);
    arb_clear(C_th);
    arb_clear(C_l);
    arb_clear(r_m);
    arb_clear(mu_th);
    arb_clear(zeta_G_r);
    arb_clear(zeta_G_r_prime);
    arb_clear(cal_c);
    
    return 0;
}

 
