#include "typical_profile.h"
#include <stdlib.h>
#include <unistd.h> //判断文件是否存在及其状态

//为了得到typical (i.e. mean) profile ζ of a given Gaußian random field ζ
//需要定义一些一系列辅助函数来计算ζ(r)


int interior_help_sigma_n(arb_t res, const arb_t k, void* params, const slong order, slong prec)
{
    //函数中所用变量
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    //梯度矩 sigma_n^2=\int dk/k * k^{2n} * P(k)
    //以 ln(k) 作为变量： sigma_n^2=\int dk * exp(k)^{2n} * P(k)
    
    // exp(k)^{2n} * P(k)
    arb_exp(s,k,2*prec); // k在30左右，提高计算精度 2*prec
    arb_pow_ui(s,s,order,2*prec); // order=2*n 直接传入
    
    power_spectrum(t,k,prec);
    arb_mul(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    
    return 0;
}


//首先是gradient moments σ_n
int Help_sigma_n_square(arb_t res, const slong n,slong prec)
{
    //函数中所用变量
    arb_t s;
    arb_init(s);
    
    slong ord;
    
    ord=n*2; // 这里传入 2*n ，避免多次重复计算
    
    // interior_help_sigma_n
    
    if(Power_spectrum_type==delta_type)
    {
        //σ_n^2=k_*^{2n}*A
        arb_pow_ui(s,K_star,ord,prec);
        arb_mul(res,s,Power_A,prec);
    }else
    {
        int ret_judge=0;
        ret_judge=Integration_arb(s, interior_help_sigma_n,NULL,ord, 
                                  Int_sigma_n_min, Int_sigma_n_max,Int_sigma_n_precision,
                                  Integration_iterate_min,Integration_iterate_max, prec);
        arb_set(res,s);
        if(ret_judge==1)
        {
            printf("Help_sigma_n_square \t %li \t 达到最大迭代次数\n", n);
        }
    }
    
    arb_clear(s);
    
    return 0;
}



//计算 sin(x)/x 及其各阶导数
int Help_sinc_n(arb_t res, const arb_t x, const slong order, slong prec)
{
    if(order==0) //原函数
    {
        arb_sinc(res,x,prec);
        
        return 0;
    }
    
    arb_t s,t,w,sin,cos;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(sin);
    arb_init(cos);
    
    arb_sin(sin,x,prec); //sin(x)
    arb_cos(cos,x,prec); //cos(x)
    
    switch(order) 
    {
        case 1 : //一阶导
            //-(sin(x)-x*cos(x))/x^2=(x*cos(x)-sin(x))/x^2
            
            arb_mul(s,x,cos,prec);
            arb_sub(s,s,sin,prec);
            arb_sqr(t,x,prec);
            arb_div(res,s,t,prec);
            
            break;
        case 2 :  //二阶导
            //-((x^2-2)*sin(x)+2*x*cos(x))/x^3
            
            arb_sqr(s,x,prec);
            arb_sub_si(s,s,2,prec);
            arb_mul(s,s,sin,prec);
            
            arb_mul_si(t,x,2,prec);
            arb_mul(t,t,cos,prec);
            arb_add(s,s,t,prec);
            
            arb_neg(s,s);
            arb_pow_ui(t,x,3,prec);
            arb_div(res,s,t,prec);
            
            
            break; 
        case 3 : //三阶导
            //((3*x^2-6)*sin(x)+(6*x-x^3)*cos(x))/x^4
            
            arb_sqr(s,x,prec);
            arb_mul_si(s,s,3,prec);
            arb_sub_si(s,s,6,prec);
            arb_mul(s,s,sin,prec);
            
            arb_mul_si(t,x,6,prec);
            arb_pow_ui(w,x,3,prec);
            arb_sub(t,t,w,prec);
            arb_mul(t,t,cos,prec);
            arb_add(s,s,t,prec);
            
            arb_pow_ui(w,x,4,prec);
            arb_div(res,s,w,prec);
            
            break;
        case 4 : //四阶导
            //((x^4-12*x^2+24)*sin(x)+(4*x^3-24*x)*cos(x))/x^5
            
            arb_pow_ui(s,x,4,prec);
            arb_sqr(t,x,prec);
            arb_mul_si(t,t,12,prec);
            arb_sub(s,s,t,prec);
            arb_add_si(s,s,24,prec);
            arb_mul(s,s,sin,prec);
            
            arb_pow_ui(t,x,3,prec);
            arb_mul_si(t,t,4,prec);
            arb_mul_si(w,x,24,prec);
            arb_sub(t,t,w,prec);
            arb_mul(t,t,cos,prec);
            arb_add(s,s,t,prec);
            
            arb_pow_ui(w,x,5,prec);
            arb_div(res,s,w,prec);
            
            break;
        default:
            printf("General -> typical_profile -> Help_sinc_n 输入有误\n");
            exit(1);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(sin);
    arb_clear(cos);
    
    return 0;
}


//计算 ψ_n(r) 及其各阶导数
int interior_help_psi_n(arb_t res, const arb_t k, void* r, const slong order, slong prec)
{
    arb_t s,t,x,exp_k;
    arb_init(s);
    arb_init(t);
    arb_init(x);
    arb_init(exp_k);
    
    //此处，半径 r 以参数作指针的形式传入
    //以 ln(k) 作为自变量
    //exp(k)^2n * sinc[exp(k)*r] * P(k)，其中n=1
    //exp(k)^2 * sin(x)/x * P(k) 其中 x=exp(k)*r ，n=1，只算 psi_1_n 
    
    arb_exp(exp_k,k,prec); //以 lnk 作变量
    arb_mul(x,exp_k,r,prec);
    
    switch(order) 
    {
        case 0 : //没有求导
            
            Help_sinc_n(s,x,0,prec); //sinc(x)
            
            break;
        case 1 : //一阶导
            
            Help_sinc_n(s,x,1,prec); //[sinc(x)]'  对r求导 x=kr
            arb_mul(s,s,exp_k,prec);
            
            break;
        case 2 : //二阶导
            Help_sinc_n(s,x,2,prec); //[sinc(x)]''  对r求导 x=kr
            arb_mul(s,s,exp_k,prec);
            arb_mul(s,s,exp_k,prec);
            
            break;
        case 3 : //三阶导
            Help_sinc_n(s,x,3,prec); //[sinc(x)]'''  对r求导 x=kr
            arb_mul(s,s,exp_k,prec);
            arb_mul(s,s,exp_k,prec);
            arb_mul(s,s,exp_k,prec);
            
            break;
        case 4 : //四阶导
            Help_sinc_n(s,x,4,prec); //[sinc(x)]''''  对r求导 x=kr
            arb_mul(s,s,exp_k,prec);
            arb_mul(s,s,exp_k,prec);
            arb_mul(s,s,exp_k,prec);
            arb_mul(s,s,exp_k,prec);
            
            break;
        default:
            printf("General -> typical_profile -> interior_help_psi_1_n 输入有误\n");
            exit(1);
    }
    
    
    if(Peak_theory_source_zeta_gradient==true)
    {
        //ψ_1(r) 前面系数 exp(k)^2
        arb_pow_ui(t,exp_k,2,prec); 
        arb_mul(s,s,t,prec);
    }else
    {
        //ψ_0 前面系数 exp(k)^0=1
        ;
    }
    
    power_spectrum(t,k,prec); //高斯情况下，ζ_G 的功率谱
    arb_mul(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(x);
    arb_clear(exp_k);
    
    return 0;
}


//计算 ψ_1(r) 及其各阶导数
int Help_psi_n(arb_t res, const arb_t r, const slong order, slong prec)
{
    //拟合暂时不用
    /*
    switch(order) 
    {
        //判断数据拟合的文件是否存在，存在直接读取相应数据，不计算
        case 0 : //没有求导
            if ( FIT_FUNC_IF && access(File_psi_1_0_fit, F_OK)==0) //存在 0，不存在 -1
            {
                if(Fit_psi_1_0==NULL)
                {
                    printf("\n发现 ψ_1(r) 拟合数据，尝试使用：\n");
                    
                    //拟合数据拟合恢复
                    Fit_psi_1_0=NULL; //利用二级指针对 Fit_psi_1_0 进行修改赋值
                    struct FUNC_FITTED_DATE** fit_temp_psi_1_0=NULL;
                    
                    fit_temp_psi_1_0=&Fit_psi_1_0;
                    Func_fit_restore(fit_temp_psi_1_0,Func_output_fitted_file,prec);
                    
                    printf("ψ_1(r) 拟合结果获取成功\n\n");
                }
                
                //由拟合结果取值
                Fit_get_value(res,r,Fit_psi_1_0,0,prec);
                
                return 0;
            }
            break;
        case 1:
            if ( FIT_FUNC_IF && access(File_psi_1_1_fit, F_OK)==0) //存在 0，不存在 -1
            {
                if(Fit_psi_1_1==NULL)
                {
                    printf("\n发现 ψ_1(r) 一阶导拟合数据，尝试使用：\n");
                    
                    //拟合数据拟合恢复
                    Fit_psi_1_1=NULL; //利用二级指针对 Fit_psi_1_1 进行修改赋值
                    struct FUNC_FITTED_DATE** fit_temp_psi_1_1=NULL;
                    
                    fit_temp_psi_1_1=&Fit_psi_1_1;
                    Func_fit_restore(fit_temp_psi_1_1,Func_output_fitted_file,prec);
                    
                    printf("ψ_1(r) 一阶导拟合结果获取成功\n\n");
                }
                
                //由拟合结果取值
                Fit_get_value(res,r,Fit_psi_1_1,0,prec);
                
                return 0;
            }
            break;
        case 2:
            if ( FIT_FUNC_IF && access(File_psi_1_2_fit, F_OK)==0) //存在 0，不存在 -1
            {
                if(Fit_psi_1_2==NULL)
                {
                    printf("\n发现 ψ_1(r) 二阶导拟合数据，尝试使用：\n");
                    
                    //拟合数据拟合恢复
                    Fit_psi_1_2=NULL; //利用二级指针对 Fit_psi_1_2 进行修改赋值
                    struct FUNC_FITTED_DATE** fit_temp_psi_1_2=NULL;
                    
                    fit_temp_psi_1_2=&Fit_psi_1_2;
                    Func_fit_restore(fit_temp_psi_1_2,Func_output_fitted_file,prec);
                    
                    printf("ψ_1(r) 二阶导拟合结果获取成功\n\n");
                }
                
                //由拟合结果取值
                Fit_get_value(res,r,Fit_psi_1_2,0,prec);
                
                return 0;
            }
            break;
        case 3:
            if ( FIT_FUNC_IF && access(File_psi_1_3_fit, F_OK)==0) //存在 0，不存在 -1
            {
                if(Fit_psi_1_3==NULL)
                {
                    printf("\n发现 ψ_1(r) 三阶导拟合数据，尝试使用：\n");
                    
                    //拟合数据拟合恢复
                    Fit_psi_1_3=NULL; //利用二级指针对 Fit_psi_1_3 进行修改赋值
                    struct FUNC_FITTED_DATE** fit_temp_psi_1_3=NULL;
                    
                    fit_temp_psi_1_3=&Fit_psi_1_3;
                    Func_fit_restore(fit_temp_psi_1_3,Func_output_fitted_file,prec);
                    
                    printf("ψ_1(r) 三阶导拟合结果获取成功\n\n");
                }
                
                //由拟合结果取值
                Fit_get_value(res,r,Fit_psi_1_3,0,prec);
                
                return 0;
            }
            break;
        case 4:
            if ( FIT_FUNC_IF && access(File_psi_1_4_fit, F_OK)==0) //存在 0，不存在 -1
            {
                if(Fit_psi_1_4==NULL)
                {
                    printf("\n发现 ψ_1(r) 四阶导拟合数据，尝试使用：\n");
                    
                    //拟合数据拟合恢复
                    Fit_psi_1_4=NULL; //利用二级指针对 Fit_psi_1_4 进行修改赋值
                    struct FUNC_FITTED_DATE** fit_temp_psi_1_4=NULL;
                    
                    fit_temp_psi_1_4=&Fit_psi_1_4;
                    Func_fit_restore(fit_temp_psi_1_4,Func_output_fitted_file,prec);
                    
                    printf("ψ_1(r) 四阶导拟合结果获取成功\n\n");
                }
                
                //由拟合结果取值
                Fit_get_value(res,r,Fit_psi_1_4,0,prec);
                return 0;
            }
            break;
        default :
            printf("General -> typical_profile ->  Help_psi_n 输入有误\n");
            exit(1);
    }
    */
    
    arb_t s,x,r_pra;
    
    arb_init(s);
    arb_init(x);
    arb_init(r_pra);
    
    //使用新的gauss_kronrod积分算法
    int ret_judge=0;
    
    
    if(Power_spectrum_type==delta_type) //δ谱下，ψ_0 = ψ_1 = sinc(x)=sinc(r*k_star)
    {
        arb_mul(x,r,K_star,prec);
        
        switch(order) 
        {
            case 0 : //没有求导
                
                Help_sinc_n(s,x,0,prec); //sinc(x)
                
                break;
            case 1 : //一阶导
                
                Help_sinc_n(s,x,1,prec); //[sinc(x)]'  对r求导 x=kr
                arb_mul(s,s,K_star,prec);
                
                break;
            case 2 : //二阶导
                Help_sinc_n(s,x,2,prec); //[sinc(x)]''  对r求导 x=kr
                arb_mul(s,s,K_star,prec);
                arb_mul(s,s,K_star,prec);
                
                break;
            case 3 : //三阶导
                Help_sinc_n(s,x,3,prec); //[sinc(x)]'''  对r求导 x=kr
                arb_mul(s,s,K_star,prec);
                arb_mul(s,s,K_star,prec);
                arb_mul(s,s,K_star,prec);
                
                break;
            case 4 : //四阶导
                Help_sinc_n(s,x,4,prec); //[sinc(x)]''''  对r求导 x=kr
                arb_mul(s,s,K_star,prec);
                arb_mul(s,s,K_star,prec);
                arb_mul(s,s,K_star,prec);
                arb_mul(s,s,K_star,prec);
                
                break;
            default:
                printf("General -> typical_profile -> interior_help_psi_1_n 输入有误\n");
                exit(1);
        } 
        
        arb_set(res,s);
        
    }else //非δ谱需要积分
    {
        arb_set(r_pra,r); //为积分传递 r 的值
        
        //计算 ψ_0 or ψ_1
        ret_judge=Integration_arb(s, interior_help_psi_n, r_pra, order,
                                  Int_sigma_n_min, Int_sigma_n_max,Int_sigma_n_precision,
                                  Integration_iterate_min,Integration_iterate_max, prec);
        
        if(Peak_theory_source_zeta_gradient==true)
        {
            //ψ_1(r)=(积分值)/(σ_1)^2
            arb_div(res,s,Sigma_1_square,prec);
        }else
        {
            //ψ_0(r)=(积分值)/(σ_0)^2
            arb_div(res,s,Sigma_0_square,prec);
        }
    }
    
    if(ret_judge==1)
    {
        printf("Help_psi_n \t %li \t 达到最大迭代次数\n", order);
    }
    
    arb_clear(s);
    arb_clear(x);
    arb_clear(r_pra);
    
    return 0;
}



//计算 ∆ψ_1(r) 及其各阶导数
int interior_help_Laplacian_psi_n(arb_t res, const arb_t k, void* r, const slong order, slong prec)
{
    arb_t s,t,x,exp_k;
    arb_init(s);
    arb_init(t);
    arb_init(x);
    arb_init(exp_k);
    
    //此处，半径 r 以参数作指针的形式传入
    //以 ln(k) 作为自变量
    //exp(k)^2n * ∆sinc[exp(k)*r] * P(k)，其中n=1，只计算 ∆ψ_1(r)
    //exp(k)^4 * [-sin(x)/x] * P(k) 其中 x=exp(k)*r ，∆sinc(x) = k^2*[-sinc(x)]
    
    
    
    arb_exp(exp_k,k,prec); //以 lnk 作变量
    arb_mul(x,exp_k,r,prec);
    
    switch(order) 
    {
        case 0 : //没有求导
            
            Help_sinc_n(s,x,0,prec); //sinc(x)
            arb_neg(s,s);
            
            break;
        case 1 : //一阶导
            
            Help_sinc_n(s,x,1,prec); //[sinc(x)]'  对r求导 x=kr
            arb_mul(s,s,exp_k,prec);
            arb_neg(s,s);
            
            break;
        case 2 : //二阶导
            Help_sinc_n(s,x,2,prec); //[sinc(x)]''  对r求导 x=kr
            arb_mul(s,s,exp_k,prec);
            arb_mul(s,s,exp_k,prec);
            arb_neg(s,s);
            
            break;
        case 3 : //三阶导
            Help_sinc_n(s,x,3,prec); //[sinc(x)]'''  对r求导 x=kr
            arb_mul(s,s,exp_k,prec);
            arb_mul(s,s,exp_k,prec);
            arb_mul(s,s,exp_k,prec);
            arb_neg(s,s);
            
            break;
        case 4 : //四阶导
            Help_sinc_n(s,x,4,prec); //[sinc(x)]''''  对r求导 x=kr
            arb_mul(s,s,exp_k,prec);
            arb_mul(s,s,exp_k,prec);
            arb_mul(s,s,exp_k,prec);
            arb_mul(s,s,exp_k,prec);
            arb_neg(s,s);
            
            break;
        default:
            printf("General -> typical_profile -> interior_Laplacian_psi_1_n 输入有误\n");
            exit(1);
    }
    
    if(Peak_theory_source_zeta_gradient==true)
    {
        //ψ_1(r) 前面系数 k^4
        arb_pow_ui(t,exp_k,4,prec);
    }else
    {
        //ψ_0(r) 前面系数 k^2
        arb_pow_ui(t,exp_k,2,prec);
    }
    
    arb_mul(s,s,t,prec);
    
    power_spectrum(t,k,prec);
    arb_mul(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(x);
    arb_clear(exp_k);
    
    return 0;
}



//计算 ∆ψ_1(r) 及其各阶导数
int Help_Laplacian_psi_n(arb_t res, const arb_t r, const slong order, slong prec)
{
    // ∆ψ_1(r) 在球坐标系中，只有径向分量 ∆ --> 2/r *d/dr + d^2/dr^2 
    // 作变量替换，x=kr, ∆ --> k^2 * (2/x *d/dx + d^2/dx^2)
    //在 ψ_1(r) 中含r只有 sinc(x) 部分，故 ∆ψ_1(r) --> ∆sinc(x)
    // ∆sinc(x) = k^2*[-sinc(x)] 
    //故后面对于 ∆ψ_1(r) 的求导全部转化为 对 sinc(x) 的求导
    
    arb_t s,r_pra;
    
    arb_init(s);
    arb_init(r_pra);
    
    arb_set(r_pra,r); //传递 r 的值
    
    int ret_judge=0;
    
    ret_judge=Integration_arb(s, interior_help_Laplacian_psi_n,r_pra, order,
                              Int_sigma_n_min, Int_sigma_n_max,Int_sigma_n_precision,
                              Integration_iterate_min,Integration_iterate_max, prec);
    
    if(Peak_theory_source_zeta_gradient==true)
    {
        //ψ_1(r)=(积分值)/(σ_1)^2
        arb_div(res,s,Sigma_1_square,prec);
    }else
    {
        //ψ_0(r)=(积分值)/(σ_0)^2
        arb_div(res,s,Sigma_0_square,prec);
    }
    
    
    if(ret_judge==1)
    {
        printf("Help_Laplacian_psi_n \t %li \t 达到最大迭代次数\n", order);
    }
    
    
    arb_clear(s);
    arb_clear(r_pra);
    
    return 0;
}



// 高斯型 ζ_G(r) 及其各阶导数
//除了μ的版本
int zeta_Gauss_profile_n_div_mu(arb_t zeta_G_r, const arb_t r, const slong order, slong prec)
{
    //函数中所用变量
    arb_t s,t,w,x;
    
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(x);
    
    arb_t phi,D_phi;//函数中所用变量
    
    arb_init(phi);
    arb_init(D_phi);
    
    switch(Power_spectrum_type) 
    {
        case lognormal_type :
        case power_law_type :
        case broken_power_law_type :
        case box_type :
        case link_cmb_type :
        case numerical_cal_type :
            
            //对应高斯型 ζ_G 的计算
            
            //包含两个参数 PT_mu、PT_k, 这里，默认 ζ_∞ = 0
            //这里的求和分成两个部分：一是 ψ_1(r)，一是 ∆ψ_1(r)
            
            if( PT_profile_simplify )
            {
                Help_psi_n(phi,r,order,prec);
                arb_set(zeta_G_r,phi);
                
                break;
            }
            
            if(Peak_theory_source_zeta_gradient==true)
            {
                // ψ_1(r) 部分
                arb_sqr(s,Gamma_3,prec); // 1-（Gamma_3）^2 
                arb_neg(s,s);
                arb_add_si(s,s,1,prec);
                
                arb_sqr(t,PT_k,prec);
                arb_mul(t,t,Gamma_3,prec);
                arb_neg(t,t);
                arb_add_si(t,t,1,prec);
                
                arb_div(t,t,s,prec); //系数完
                Help_psi_n(phi,r,order,prec);
                arb_mul(phi,t,phi,prec);
                
                // ∆ψ_1(r) 部分
                arb_sqr(t,PT_k,prec);
                arb_div(t,t,Gamma_3,prec);
                arb_neg(t,t);
                arb_add_si(t,t,1,prec);
                arb_div(t,t,s,prec); //s用掉
                
                arb_sqr(s,R_3,prec);
                arb_div_si(s,s,3,prec);
                arb_mul(t,t,s,prec); //系数完
                
                Help_Laplacian_psi_n(D_phi,r,order,prec);
                arb_mul(D_phi,t,D_phi,prec);
                
                arb_add(zeta_G_r,phi,D_phi,prec); //右边完
                
            }else
            {
                //这里作了变量替换 \tilde{k}_1 ^ = k_1^2 * σ_0/σ_2
                // ψ_0(r) 部分
                arb_sqr(s,Gamma_1,prec); // 1-（Gamma_1）^2 
                arb_neg(s,s);
                arb_add_si(s,s,1,prec);
                
                arb_sqr(t,PT_k,prec);
                arb_mul(t,t,Gamma_1,prec);
                arb_neg(t,t);
                arb_add_si(t,t,1,prec);
                
                arb_div(t,t,s,prec); //系数完
                Help_psi_n(phi,r,order,prec);
                arb_mul(phi,t,phi,prec);
                
                // ∆ψ_0(r) 部分
                arb_sqr(t,PT_k,prec);
                arb_div(t,t,Gamma_1,prec);
                arb_neg(t,t);
                arb_add_si(t,t,1,prec);
                arb_div(t,t,s,prec); //s用掉
                
                arb_sqr(s,R_1,prec);
                arb_div_si(s,s,3,prec);
                arb_mul(t,t,s,prec); //系数完
                
                Help_Laplacian_psi_n(D_phi,r,order,prec);
                arb_mul(D_phi,t,D_phi,prec);
                
                arb_add(zeta_G_r,phi,D_phi,prec); //右边完
                
            }
            
            
            
            break;
            
            case delta_type :
                
                //对应高斯型 ζ_G 的计算
                //注意到，对于δ的情况，无论是否采用梯度，得到的结果一样
                
                //首先计算 x
                arb_mul(x,K_star,r,prec); //x=k_star*r
                
                switch(order)
                {
                    case 0: //原函数
                        //功率谱为delta函数时，ζ(r) 易解析求出
                        // ζ(r)=μ*sinc(k_star*r)
                        
                        Help_sinc_n(s,x,0,prec);
                        arb_set(zeta_G_r,s);
                        break;
                        
                    case 1: //一阶导
                        //功率谱为delta函数时易解析求出 x=k_star*r
                        //ζ(r)=μ*sinc(x)
                        //ζ(r)^prime=μ * sinc'(x) * k_star
                        
                        Help_sinc_n(s,x,1,prec);
                        arb_mul(zeta_G_r,s,K_star,prec); //一阶导掉个 k_star
                        
                        break;
                        
                    case 2: //二阶导
                        //功率谱为delta函数时易解析求出 x=k_star*r
                        //ζ(r)=μ*sinc(x)
                        //ζ(r)^prime=μ * 1/x *[cos(x)-sinc(x)] * k_star
                        //ζ(r)^prime^prime=μ * (-1/x^2) * [x*sin(x)+2(cos(x)-sinc(x))] * (k_star)^2
                        
                        Help_sinc_n(s,x,2,prec);
                        arb_mul(s,s,K_star,prec); //二阶导掉个 (k_star)^2
                        arb_mul(zeta_G_r,s,K_star,prec);
                        
                        break;
                        
                    case 3: //三阶导
                        //ζ(r)=μ*sinc(x)  x=k_star*r
                        //求 ζ'''(r)
                        
                        Help_sinc_n(s,x,3,prec);
                        arb_mul(s,s,K_star,prec); //三阶导掉个 (k_star)^3
                        arb_mul(s,s,K_star,prec);
                        arb_mul(zeta_G_r,s,K_star,prec);
                        
                        break;
                    case 4: //四阶导
                        //ζ(r)=μ*sinc(x)   x=k_star*r
                        //求 ζ''''(r)
                        
                        Help_sinc_n(s,x,4,prec);
                        arb_mul(s,s,K_star,prec); //三阶导掉个 (k_star)^4
                        arb_mul(s,s,K_star,prec);
                        arb_mul(s,s,K_star,prec);
                        arb_mul(zeta_G_r,s,K_star,prec);
                        
                        break;
                    default:
                        printf("General->typical_profile->zeta_Gauss_profile_n_div_mu->delta_type->order 输入有误\n");
                        exit(1);
                }
                
                
                break;
                
                    default :
                        printf("General-> typical_profile->zeta_Gauss_profile_n_div_mu 中的 Power_spectrum_type 有误\n");
                        exit(1);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(x);
    
    arb_clear(phi);
    arb_clear(D_phi);
    
    return 0;
}

//没除μ的完整版
int zeta_Gauss_profile_n(arb_t res, const arb_t r, const slong order, slong prec)
{
    arb_t zeta_G_r;
    arb_init(zeta_G_r);
    
    zeta_Gauss_profile_n_div_mu(zeta_G_r,r,order,prec);
    arb_mul(res,zeta_G_r,PT_mu,prec); //ζ(r) 完
    
    arb_clear(zeta_G_r);
    
    return 0;
}



//利用 ζ=F(ζ_G) 得到的 profile ζ(r) 及其各阶导数
static void zeta_profile_n_Func_direct(arb_t res, const arb_t r, const slong order, slong prec)
{
    //函数中所用变量
    arb_t s,t,w;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    arb_t G_0,G_1,G_2,G_3;
    arb_init(G_0);
    arb_init(G_1);
    arb_init(G_2);
    arb_init(G_3);
    
    // γ_n 的大小在 1 左右，基本不受 K_star 影响
    //当 k_square等于其平均值时， ζ(r)=PT_mu * ψ_n(r)
    // 而ζ(r)中的 ∆ψ_1(r) 部分，被 (R_3)^2 极大的压低，当 K_star 非常大时，基本可忽略
    
    //这里，不用对功率谱分类，统一对 Zeta_type 分类即可
    
    switch(Zeta_type) 
    {
        case gaussian_type :
            //直接调用纯高斯情况函数
            
            zeta_Gauss_profile_n(res,r,order,prec);
            
            break;
            
        case exponential_tail_type :
            // ζ=f(ζ_G) 为一复合函数，当有导数时，需按复合函数求导规则
            switch(order)
            {
                case 0 : //原函数
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    
                    Non_Gaussianity_exponential_tail_n(res,G_0,0,prec);
                    
                    break;
                case 1 : //一阶导
                    //ζ'=f'(ζ_G)*ζ'_G
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    zeta_Gauss_profile_n(G_1,r,1,prec);
                    
                    Non_Gaussianity_exponential_tail_n(s,G_0,1,prec);
                    
                    arb_mul(res,s,G_1,prec);
                    
                    break;
                case 2 : //二阶导
                    //ζ''=f''(ζ_G)*(ζ'_G)^2 + f'(ζ_G)*ζ''_G
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    zeta_Gauss_profile_n(G_1,r,1,prec);
                    zeta_Gauss_profile_n(G_2,r,2,prec);
                    
                    //前半部分
                    
                    Non_Gaussianity_exponential_tail_n(s,G_0,2,prec);
                    arb_sqr(t,G_1,prec);
                    arb_mul(s,s,t,prec);
                    
                    //后半部分
                    Non_Gaussianity_exponential_tail_n(t,G_0,1,prec);
                    arb_mul(t,t,G_2,prec);
                    
                    arb_add(res,s,t,prec);
                    
                    break;
                case 3 : //三阶导
                    //ζ''=f''(ζ_G)*(ζ'_G)^2 + f'(ζ_G)*ζ''_G
                    //ζ''' = f'''(ζ_G)*(ζ'_G)^3 + f''(ζ_G)*2*(ζ'_G)*(ζ''_G)
                    //       + f''(ζ_G)*(ζ'_G)*ζ''_G + f'(ζ_G)*ζ'''_G
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    zeta_Gauss_profile_n(G_1,r,1,prec);
                    zeta_Gauss_profile_n(G_2,r,2,prec);
                    zeta_Gauss_profile_n(G_3,r,3,prec);
                    
                    //前半部分
                    Non_Gaussianity_exponential_tail_n(s,G_0,3,prec);
                    arb_pow_ui(t,G_1,3,prec);
                    arb_mul(s,s,t,prec);
                    
                    //中间
                    Non_Gaussianity_exponential_tail_n(t,G_0,2,prec);
                    arb_mul(w,G_1,G_2,prec);
                    arb_mul_si(w,w,3,prec); //中间两项可以合并
                    arb_mul(t,t,w,prec);
                    arb_add(s,s,t,prec);
                    
                    //后半部分
                    Non_Gaussianity_exponential_tail_n(t,G_0,1,prec);
                    arb_mul(t,t,G_3,prec);
                    
                    arb_add(res,s,t,prec);
                    
                    break;
                case 4 : //四阶导
                    printf("General -> typical_profile -> zeta_profile_n->lognormal_type->Zeta_type->exponential_tail_type 阶数n输入有误\n");
                    exit(1);
                    break;
                default :
                    printf("General -> typical_profile -> zeta_profile_n->lognormal_type->Zeta_type->exponential_tail_type 阶数n输入有误\n");
                    exit(1);
            }
            
            break;
            
        case up_step_type :
            // ζ=f(ζ_G) 为一复合函数，当有导数时，需按复合函数求导规则
            switch(order)
            {
                case 0 : //原函数
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    
                    Non_Gaussianity_up_step_n(res,G_0,0,prec);
                    
                    break;
                case 1 : //一阶导
                    //ζ'=f'(ζ_G)*ζ'_G
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    zeta_Gauss_profile_n(G_1,r,1,prec);
                    
                    Non_Gaussianity_up_step_n(s,G_0,1,prec);
                    
                    arb_mul(res,s,G_1,prec);
                    break;
                case 2 : //二阶导
                    //ζ''=f''(ζ_G)*(ζ'_G)^2 + f'(ζ_G)*ζ''_G
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    zeta_Gauss_profile_n(G_1,r,1,prec);
                    zeta_Gauss_profile_n(G_2,r,2,prec);
                    
                    //前半部分
                    
                    Non_Gaussianity_up_step_n(s,G_0,2,prec);
                    arb_sqr(t,G_1,prec);
                    arb_mul(s,s,t,prec);
                    
                    //后半部分
                    Non_Gaussianity_up_step_n(t,G_0,1,prec);
                    arb_mul(t,t,G_2,prec);
                    
                    arb_add(res,s,t,prec);
                    break;
                case 3 : //三阶导
                    //ζ''=f''(ζ_G)*(ζ'_G)^2 + f'(ζ_G)*ζ''_G
                    //ζ''' = f'''(ζ_G)*(ζ'_G)^3 + f''(ζ_G)*2*(ζ'_G)*(ζ''_G)
                    //       + f''(ζ_G)*(ζ'_G)*ζ''_G + f'(ζ_G)*ζ'''_G
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    zeta_Gauss_profile_n(G_1,r,1,prec);
                    zeta_Gauss_profile_n(G_2,r,2,prec);
                    zeta_Gauss_profile_n(G_3,r,3,prec);
                    
                    //前半部分
                    Non_Gaussianity_up_step_n(s,G_0,3,prec);
                    arb_pow_ui(t,G_1,3,prec);
                    arb_mul(s,s,t,prec);
                    
                    //中间
                    Non_Gaussianity_up_step_n(t,G_0,2,prec);
                    arb_mul(w,G_1,G_2,prec);
                    arb_mul_si(w,w,3,prec); //中间两项可以合并
                    arb_mul(t,t,w,prec);
                    arb_add(s,s,t,prec);
                    
                    //后半部分
                    Non_Gaussianity_up_step_n(t,G_0,1,prec);
                    arb_mul(t,t,G_3,prec);
                    
                    arb_add(res,s,t,prec);
                    break;
                case 4 : //四阶导
                    printf("General -> typical_profile -> zeta_profile_n-->Zeta_type->up_step 阶数n输入有误\n");
                    exit(1);
                    break;
                default :
                    printf("General -> typical_profile -> zeta_profile_n-->Zeta_type->up_step 阶数n输入有误\n");
                    exit(1);
            }
            
            break;
            
        case power_expansion_type :
            // ζ=f(ζ_G) 为一复合函数，当有导数时，需按复合函数求导规则
            switch(order)
            {
                case 0 : //原函数
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    
                    Non_Gaussianity_power_expansion_n(res,G_0,0,prec);
                    
                    break;
                case 1 : //一阶导
                    //ζ'=f'(ζ_G)*ζ'_G
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    zeta_Gauss_profile_n(G_1,r,1,prec);
                    
                    Non_Gaussianity_power_expansion_n(s,G_0,1,prec);
                    
                    arb_mul(res,s,G_1,prec);
                    break;
                case 2 : //二阶导
                    //ζ''=f''(ζ_G)*(ζ'_G)^2 + f'(ζ_G)*ζ''_G
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    zeta_Gauss_profile_n(G_1,r,1,prec);
                    zeta_Gauss_profile_n(G_2,r,2,prec);
                    
                    //前半部分
                    
                    Non_Gaussianity_power_expansion_n(s,G_0,2,prec);
                    arb_sqr(t,G_1,prec);
                    arb_mul(s,s,t,prec);
                    
                    //后半部分
                    Non_Gaussianity_power_expansion_n(t,G_0,1,prec);
                    arb_mul(t,t,G_2,prec);
                    
                    arb_add(res,s,t,prec);
                    break;
                case 3 : //三阶导
                    //ζ''=f''(ζ_G)*(ζ'_G)^2 + f'(ζ_G)*ζ''_G
                    //ζ''' = f'''(ζ_G)*(ζ'_G)^3 + f''(ζ_G)*2*(ζ'_G)*(ζ''_G)
                    //       + f''(ζ_G)*(ζ'_G)*ζ''_G + f'(ζ_G)*ζ'''_G
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    zeta_Gauss_profile_n(G_1,r,1,prec);
                    zeta_Gauss_profile_n(G_2,r,2,prec);
                    zeta_Gauss_profile_n(G_3,r,3,prec);
                    
                    //前半部分
                    Non_Gaussianity_power_expansion_n(s,G_0,3,prec);
                    arb_pow_ui(t,G_1,3,prec);
                    arb_mul(s,s,t,prec);
                    
                    //中间
                    Non_Gaussianity_power_expansion_n(t,G_0,2,prec);
                    arb_mul(w,G_1,G_2,prec);
                    arb_mul_si(w,w,3,prec); //中间两项可以合并
                    arb_mul(t,t,w,prec);
                    arb_add(s,s,t,prec);
                    
                    //后半部分
                    Non_Gaussianity_power_expansion_n(t,G_0,1,prec);
                    arb_mul(t,t,G_3,prec);
                    
                    arb_add(res,s,t,prec);
                    break;
                case 4 : //四阶导
                    printf("General -> typical_profile -> zeta_profile_n-->Zeta_type-->power_expansion 阶数n输入有误\n");
                    exit(1);
                    break;
                default :
                    printf("General -> typical_profile -> zeta_profile_n-->Zeta_type->power_expansion 阶数n输入有误\n");
                    exit(1);
            }
            
            break;
        case narrow_step_1_type :
            // ζ=f(ζ_G) 为一复合函数，当有导数时，需按复合函数求导规则
            switch(order)
            {
                case 0 : //原函数
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    
                    Non_Gaussianity_narrow_1_up_step_n(res,G_0,0,prec);
                    break;
                case 1 : //一阶导
                    //ζ'=f'(ζ_G)*ζ'_G
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    zeta_Gauss_profile_n(G_1,r,1,prec);
                    
                    Non_Gaussianity_narrow_1_up_step_n(s,G_0,1,prec);
                    
                    arb_mul(res,s,G_1,prec);
                    break;
                case 2 : //二阶导
                    //ζ''=f''(ζ_G)*(ζ'_G)^2 + f'(ζ_G)*ζ''_G
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    zeta_Gauss_profile_n(G_1,r,1,prec);
                    zeta_Gauss_profile_n(G_2,r,2,prec);
                    
                    //前半部分
                    
                    Non_Gaussianity_narrow_1_up_step_n(s,G_0,2,prec);
                    arb_sqr(t,G_1,prec);
                    arb_mul(s,s,t,prec);
                    
                    //后半部分
                    Non_Gaussianity_narrow_1_up_step_n(t,G_0,1,prec);
                    arb_mul(t,t,G_2,prec);
                    
                    arb_add(res,s,t,prec);
                    break;
                case 3 : //三阶导
                    //ζ''=f''(ζ_G)*(ζ'_G)^2 + f'(ζ_G)*ζ''_G
                    //ζ''' = f'''(ζ_G)*(ζ'_G)^3 + f''(ζ_G)*2*(ζ'_G)*(ζ''_G)
                    //       + f''(ζ_G)*(ζ'_G)*ζ''_G + f'(ζ_G)*ζ'''_G
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    zeta_Gauss_profile_n(G_1,r,1,prec);
                    zeta_Gauss_profile_n(G_2,r,2,prec);
                    zeta_Gauss_profile_n(G_3,r,3,prec);
                    
                    //前半部分
                    Non_Gaussianity_narrow_1_up_step_n(s,G_0,3,prec);
                    arb_pow_ui(t,G_1,3,prec);
                    arb_mul(s,s,t,prec);
                    
                    //中间
                    Non_Gaussianity_narrow_1_up_step_n(t,G_0,2,prec);
                    arb_mul(w,G_1,G_2,prec);
                    arb_mul_si(w,w,3,prec); //中间两项可以合并
                    arb_mul(t,t,w,prec);
                    arb_add(s,s,t,prec);
                    
                    //后半部分
                    Non_Gaussianity_narrow_1_up_step_n(t,G_0,1,prec);
                    arb_mul(t,t,G_3,prec);
                    
                    arb_add(res,s,t,prec);
                    break;
                case 4 : //四阶导
                    printf("General -> typical_profile -> zeta_profile_n-->Zeta_type->narrow_step_1_type 阶数n输入有误\n");
                    exit(1);
                    break;
                default :
                    printf("General -> typical_profile -> zeta_profile_n-->Zeta_type->narrow_step_1_type 阶数n输入有误\n");
                    exit(1);
            }
            
            break;
         case narrow_step_1_2_type :
            // ζ=f(ζ_G) 为一复合函数，当有导数时，需按复合函数求导规则
            switch(order)
            {
                case 0 : //原函数
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    
                    Non_Gaussianity_narrow_1_2_up_step_n(res,G_0,0,prec);
                    
                    break;
                case 1 : //一阶导
                    //ζ'=f'(ζ_G)*ζ'_G
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    zeta_Gauss_profile_n(G_1,r,1,prec);
                    
                    Non_Gaussianity_narrow_1_2_up_step_n(s,G_0,1,prec);
                    
                    arb_mul(res,s,G_1,prec);
                    break;
                case 2 : //二阶导
                    //ζ''=f''(ζ_G)*(ζ'_G)^2 + f'(ζ_G)*ζ''_G
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    zeta_Gauss_profile_n(G_1,r,1,prec);
                    zeta_Gauss_profile_n(G_2,r,2,prec);
                    
                    //前半部分
                    
                    Non_Gaussianity_narrow_1_2_up_step_n(s,G_0,2,prec);
                    arb_sqr(t,G_1,prec);
                    arb_mul(s,s,t,prec);
                    
                    //后半部分
                    Non_Gaussianity_narrow_1_2_up_step_n(t,G_0,1,prec);
                    arb_mul(t,t,G_2,prec);
                    
                    arb_add(res,s,t,prec);
                    break;
                case 3 : //三阶导
                    //ζ''=f''(ζ_G)*(ζ'_G)^2 + f'(ζ_G)*ζ''_G
                    //ζ''' = f'''(ζ_G)*(ζ'_G)^3 + f''(ζ_G)*2*(ζ'_G)*(ζ''_G)
                    //       + f''(ζ_G)*(ζ'_G)*ζ''_G + f'(ζ_G)*ζ'''_G
                    
                    zeta_Gauss_profile_n(G_0,r,0,prec);
                    zeta_Gauss_profile_n(G_1,r,1,prec);
                    zeta_Gauss_profile_n(G_2,r,2,prec);
                    zeta_Gauss_profile_n(G_3,r,3,prec);
                    
                    //前半部分
                    Non_Gaussianity_narrow_1_2_up_step_n(s,G_0,3,prec);
                    arb_pow_ui(t,G_1,3,prec);
                    arb_mul(s,s,t,prec);
                    
                    //中间
                    Non_Gaussianity_narrow_1_2_up_step_n(t,G_0,2,prec);
                    arb_mul(w,G_1,G_2,prec);
                    arb_mul_si(w,w,3,prec); //中间两项可以合并
                    arb_mul(t,t,w,prec);
                    arb_add(s,s,t,prec);
                    
                    //后半部分
                    Non_Gaussianity_narrow_1_2_up_step_n(t,G_0,1,prec);
                    arb_mul(t,t,G_3,prec);
                    
                    arb_add(res,s,t,prec);
                    break;
                case 4 : //四阶导
                    printf("General -> typical_profile -> zeta_profile_n-->Zeta_type->narrow_step_1_2_type 阶数n输入有误\n");
                    exit(1);
                    break;
                default :
                    printf("General -> typical_profile -> zeta_profile_n-->Zeta_type->narrow_step_1_2_type 阶数n输入有误\n");
                    exit(1);
            }
            
            break;
        default:
            printf("General -> typical_profile -> zeta_profile_n 中 zeta_type 输入有误\n");
            exit(1);
    }
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    
    arb_clear(G_0);
    arb_clear(G_1);
    arb_clear(G_2);
    arb_clear(G_3);
    
}


//
//注意到，对于非高斯性，常见的 typical profile ζ(r)=F[ζ_G(r)] 存在问题
//可以利用幂级数展开，通过高斯的功率谱 P_ζ_G(k)，计算得到非高斯的功率谱 P_ζ(k)
//再利用 P_ζ(k)，参照高斯 profile ζ_G(r) 的求法，得到非高斯的 profile ζ(r)
//


//计算非高斯修正 ψ_n(r) 及其各阶导数
//这里的函数，并不能与前面共用，后面还会有调用前面的函数，这里单独写
int interior_help_psi_n_Non_Gaussian_correction(arb_t res, const arb_t k, void* r, const slong order, slong prec)
{
    arb_t s,t,x,exp_k;
    arb_init(s);
    arb_init(t);
    arb_init(x);
    arb_init(exp_k);
    
    //此处，半径 r 以参数作指针的形式传入
    //以 ln(k) 作为自变量
    //exp(k)^2n * sinc[exp(k)*r] * P(k)
    //如 n=1 时，exp(k)^2 * sin(x)/x * P(k) 其中 x=exp(k)*r 
    
    arb_exp(exp_k,k,prec); //以 lnk 作变量
    arb_mul(x,exp_k,r,prec);
    
    switch(order) 
    {
        case 0 : //没有求导
            
            Help_sinc_n(s,x,0,prec); //sinc(x)
            
            break;
        case 1 : //一阶导
            
            Help_sinc_n(s,x,1,prec); //[sinc(x)]'  对r求导 x=kr
            arb_mul(s,s,exp_k,prec);
            
            break;
        case 2 : //二阶导
            Help_sinc_n(s,x,2,prec); //[sinc(x)]''  对r求导 x=kr
            arb_mul(s,s,exp_k,prec);
            arb_mul(s,s,exp_k,prec);
            
            break;
        case 3 : //三阶导
            Help_sinc_n(s,x,3,prec); //[sinc(x)]'''  对r求导 x=kr
            arb_mul(s,s,exp_k,prec);
            arb_mul(s,s,exp_k,prec);
            arb_mul(s,s,exp_k,prec);
            
            break;
        case 4 : //四阶导
            Help_sinc_n(s,x,4,prec); //[sinc(x)]''''  对r求导 x=kr
            arb_mul(s,s,exp_k,prec);
            arb_mul(s,s,exp_k,prec);
            arb_mul(s,s,exp_k,prec);
            arb_mul(s,s,exp_k,prec);
            
            break;
        default:
            printf("General -> typical_profile -> interior_help_psi_n_Non_Gaussian_correction 输入有误\n");
            exit(1);
    }
    
    //我们这里实际计算的是，n=0
    //ψ_0 前面系数 exp(k)^0=1，故这里可以什么也不干
    
    power_spectrum_non_Gaussian(t,k,prec); //非高斯情况下，ζ 的功率谱
    arb_mul(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(x);
    arb_clear(exp_k);
    
    return 0;
}


//计算非高斯修正 ψ_1(r) 及其各阶导数
//这里的函数，并不能与前面共用，后面还会有调用前面的函数，这里单独写
int Help_psi_n_Non_Gaussian_correction(arb_t res, const arb_t r, const slong order, slong prec)
{
    arb_t s,w,x,r_pra;
    
    arb_init(s);
    arb_init(w);
    arb_init(x);
    arb_init(r_pra);
    
    int ret_judge=0;
    
    if(Power_spectrum_type==delta_type) //δ谱下，ψ_0 = ψ_1 = sinc(x)=sinc(r*k_star)
    {
        arb_mul(x,r,K_star,prec);
        
        switch(order) 
        {
            case 0 : //没有求导
                
                Help_sinc_n(s,x,0,prec); //sinc(x)
                
                break;
            case 1 : //一阶导
                
                Help_sinc_n(s,x,1,prec); //[sinc(x)]'  对r求导 x=kr
                arb_mul(s,s,K_star,prec);
                
                break;
            case 2 : //二阶导
                Help_sinc_n(s,x,2,prec); //[sinc(x)]''  对r求导 x=kr
                arb_mul(s,s,K_star,prec);
                arb_mul(s,s,K_star,prec);
                
                break;
            case 3 : //三阶导
                Help_sinc_n(s,x,3,prec); //[sinc(x)]'''  对r求导 x=kr
                arb_mul(s,s,K_star,prec);
                arb_mul(s,s,K_star,prec);
                arb_mul(s,s,K_star,prec);
                
                break;
            case 4 : //四阶导
                Help_sinc_n(s,x,4,prec); //[sinc(x)]''''  对r求导 x=kr
                arb_mul(s,s,K_star,prec);
                arb_mul(s,s,K_star,prec);
                arb_mul(s,s,K_star,prec);
                arb_mul(s,s,K_star,prec);
                
                break;
            default:
                printf("General -> typical_profile -> Help_psi_n_Non_Gaussian_correction 输入有误\n");
                exit(1);
        } 
        
        arb_set(res,s);
        
    }else //非δ谱需要积分
    {
        arb_set(r_pra,r); //为积分传递 r 的值
        
        //这里的积分精度，需根据阶数动态调整，每求一次会掉个 k_star
        //而这里的精度是绝对精度
        
        if(order==0)
        {
            arb_set_str(w,"1E-10",prec);
        }else if(order==1)
        {
           arb_set_str(w,"1E-7",prec);
        }else if (order==2)
        {
            arb_set_str(w,"1E-2",prec);
        }else if (order==3)
        {
            arb_set_str(w,"1",prec);
        }else if (order==4)
        {
            arb_set_str(w,"1E3",prec);
        }
        
        //计算 ψ_0 or ψ_1
        ret_judge=Integration_arb(s, interior_help_psi_n_Non_Gaussian_correction, r_pra, order,
                                  Int_sigma_n_min, Int_sigma_n_max,w,
                                  4,8, prec);
        
        //我们这里实际计算的是，n=0
        //ψ_0(r)=(积分值)/(σ_0)^2
        arb_div(res,s,Sigma_0_square,prec);
    }
    
    if(ret_judge==1)
    {
        //printf("Help_psi_n_Non_Gaussian_correction \t %li \t 达到最大迭代次数\n", order);
    }
    
    arb_clear(s);
    arb_clear(w);
    arb_clear(x);
    arb_clear(r_pra);
    
    return 0;
}


//利用 Non_Gaussian 功率谱 P_ζ(k) 得到的 profile ζ(r) 及其各阶导数
void zeta_profile_n_Non_Gaussian_correction(arb_t res, const arb_t r, const slong order, slong prec)
{
    //函数中所用变量
    arb_t s,t,x;
    
    arb_init(s);
    arb_init(t);
    arb_init(x);
    
    arb_t phi,zeta_G_r;
    
    arb_init(phi);
    arb_init(zeta_G_r);
    
    
    switch(Power_spectrum_type) 
    {
        case lognormal_type :
        case power_law_type :
        case broken_power_law_type :
        case box_type :
        case link_cmb_type :
        case numerical_cal_type :
            
            //对应高斯型 ζ_G 的计算
            Help_psi_n_Non_Gaussian_correction(phi,r,order,prec);
            arb_set(zeta_G_r,phi);
            
            break;
        case delta_type :
            
            //对应高斯型 ζ_G 的计算
            //注意到，对于δ的情况，无论是否采用梯度，得到的结果一样
            
            //首先计算 x
            arb_mul(x,K_star,r,prec); //x=k_star*r
            
            switch(order)
            {
                case 0: //原函数
                    //功率谱为delta函数时，ζ(r) 易解析求出
                    // ζ(r)=μ*sinc(k_star*r)
                    
                    Help_sinc_n(s,x,0,prec);
                    arb_set(zeta_G_r,s);
                    break;
                    
                case 1: //一阶导
                    //功率谱为delta函数时易解析求出 x=k_star*r
                    //ζ(r)=μ*sinc(x)
                    //ζ(r)^prime=μ * sinc'(x) * k_star
                    
                    Help_sinc_n(s,x,1,prec);
                    arb_mul(zeta_G_r,s,K_star,prec); //一阶导掉个 k_star
                    
                    break;
                    
                case 2: //二阶导
                    //功率谱为delta函数时易解析求出 x=k_star*r
                    //ζ(r)=μ*sinc(x)
                    //ζ(r)^prime=μ * 1/x *[cos(x)-sinc(x)] * k_star
                    //ζ(r)^prime^prime=μ * (-1/x^2) * [x*sin(x)+2(cos(x)-sinc(x))] * (k_star)^2
                    
                    Help_sinc_n(s,x,2,prec);
                    arb_mul(s,s,K_star,prec); //二阶导掉个 (k_star)^2
                    arb_mul(zeta_G_r,s,K_star,prec);
                    
                    break;
                    
                case 3: //三阶导
                    //ζ(r)=μ*sinc(x)  x=k_star*r
                    //求 ζ'''(r)
                    
                    Help_sinc_n(s,x,3,prec);
                    arb_mul(s,s,K_star,prec); //三阶导掉个 (k_star)^3
                    arb_mul(s,s,K_star,prec);
                    arb_mul(zeta_G_r,s,K_star,prec);
                    
                    break;
                case 4: //四阶导
                    //ζ(r)=μ*sinc(x)   x=k_star*r
                    //求 ζ''''(r)
                    
                    Help_sinc_n(s,x,4,prec);
                    arb_mul(s,s,K_star,prec); //三阶导掉个 (k_star)^4
                    arb_mul(s,s,K_star,prec);
                    arb_mul(s,s,K_star,prec);
                    arb_mul(zeta_G_r,s,K_star,prec);
                    
                    break;
                default:
                    printf("General->typical_profile->zeta_profile_n_Non_Gaussian_correction->delta_type->order 输入有误\n");
                    exit(1);
            }
            
            break;
        default :
                printf("General-> typical_profile->zeta_profile_n_Non_Gaussian_correction 中的  Power_spectrum_type 有误\n");
                exit(1);
    }
    
    
    arb_mul(res,zeta_G_r,PT_mu,prec); //ζ(r) 完
    
    switch(Zeta_type) 
    {
        case up_step_type :
            //这里还需要对 ζ(r) 的大小作判断，ζ < 2/|h|
            
            arb_abs(s,Up_step_h); // 2/|h|
            arb_inv(s,s,prec);
            arb_mul_ui(s,s,2,prec);
            
            Help_psi_n_Non_Gaussian_correction(t,r,0,prec); //ζ(r)
            arb_mul(t,t,PT_mu,prec);
            
            if( arb_gt(t,s) )
            {
                arb_zero(res);
            }
            
            break;
        default:
            printf("General -> typical_profile -> zeta_profile_n_Non_Gaussian_correction 中 zeta_type 输入有误\n");
            exit(1);
    }
    
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(x);
    
    arb_clear(phi);
    arb_clear(zeta_G_r);
}


int zeta_profile_n(arb_t res, const arb_t r, const slong order, slong prec) // ζ(r) 及其各阶导数
{
    
    if( Non_Gaussian_typical_profile_correction==true )
    {
        zeta_profile_n_Non_Gaussian_correction(res, r, order, prec); //利用NG功率谱P_ζ(k)得到的profile ζ(r)及其各阶导数
    }else
    {
        zeta_profile_n_Func_direct(res, r, order, prec); //利用 ζ=F(ζ_G) 得到的 profile ζ(r) 及其各阶导数
    }
    
    return 0;
}
