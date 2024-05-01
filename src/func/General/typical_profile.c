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
    slong ord;
    
    ord=n*2; // 这里传入 2*n ，避免多次重复计算
    
    // interior_help_sigma_n
    
    int ret_judge=0;
    ret_judge=Integration_arb(res, interior_help_sigma_n,NULL,ord, 
                                        Int_sigma_n_min, Int_sigma_n_max,Int_sigma_n_precision,
                                        Integration_iterate_min,Integration_iterate_max, prec);
    if(ret_judge==1)
    {
        printf("Help_sigma_n_square \t %li \t 达到最大迭代次数\n", n);
    }
    
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


//计算 ψ_1(r) 及其各阶导数
int interior_help_psi_1_n(arb_t res, const arb_t k, void* r, const slong order, slong prec)
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
    
    arb_pow_ui(t,exp_k,2,prec); //前面系数 exp(k)^2 ，只算 psi_1_n 
    arb_mul(s,s,t,prec);
    
    power_spectrum(t,k,prec);
    arb_mul(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(x);
    arb_clear(exp_k);
    
    return 0;
}

//计算 ψ_1(r) 及其各阶导数
int Help_psi_1_n(arb_t res, const arb_t r, const slong order, slong prec)
{
    
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
            printf("General -> typical_profile ->  Help_psi_1_n 输入有误\n");
            exit(1);
    }
    
    
    arb_t s,r_pra;
    
    arb_init(s);
    arb_init(r_pra);
    
    arb_set(r_pra,r); //传递 r 的值
    
    //使用新的gauss_kronrod积分算法
    int ret_judge=0;
    ret_judge=Integration_arb(s, interior_help_psi_1_n, r_pra, order,
                                        Int_sigma_n_min, Int_sigma_n_max,Int_sigma_n_precision,
                                        Integration_iterate_min,Integration_iterate_max, prec);
    if(ret_judge==1)
    {
        printf("Help_psi_1_n \t %li \t 达到最大迭代次数\n", order);
    }
    
    
    //ψ_1(r)=(积分值)/(σ_1)^2
    arb_div(res,s,Sigma_1_square,prec);
    
    arb_clear(s);
    arb_clear(r_pra);
    
    return 0;
}


//计算 ∆ψ_1(r) 及其各阶导数
int interior_help_Laplacian_psi_1_n(arb_t res, const arb_t k, void* r, const slong order, slong prec)
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
    
    arb_pow_ui(t,exp_k,4,prec); //前面系数
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
int Help_Laplacian_psi_1_n(arb_t res, const arb_t r, const slong order, slong prec)
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
    ret_judge=Integration_arb(s, interior_help_Laplacian_psi_1_n,r_pra, order,
                                        Int_sigma_n_min, Int_sigma_n_max,Int_sigma_n_precision,
                                        Integration_iterate_min,Integration_iterate_max, prec);
    if(ret_judge==1)
    {
        printf("Help_Laplacian_psi_1_n \t %li \t 达到最大迭代次数\n", order);
    }
    
    //ψ_1(r)=(积分值)/(σ_1)^2
    arb_div(res,s,Sigma_1_square,prec);
    
    arb_clear(s);
    arb_clear(r_pra);
    
    return 0;
}




// 高斯型 ζ_G(r) 及其各阶导数
int zeta_Gauss_profile_n(arb_t res, const arb_t r, const slong order, slong prec)
{
    //函数中所用变量
    arb_t s,t,x;
    
    arb_init(s);
    arb_init(t);
    arb_init(x);
     
    arb_t phi_1,D_phi_1;//函数中所用变量
    
    arb_init(phi_1);
    arb_init(D_phi_1);
    
    switch(Power_spectrum_type) 
    {
        case lognormal_type :
            
            //对应高斯型 ζ_G 的计算
            
            //包含两个参数 Mu_2、K_3_square, 这里，默认 ζ_∞ = 0
            //这里的求和分成两个部分：一是 ψ_1(r)，一是 ∆ψ_1(r)
            
            if( arb_equal(K_3_square,Gamma_3) )
            {
                Help_psi_1_n(phi_1,r,order,prec);
                
                arb_mul(res,phi_1,Mu_2,prec); // ζ(r) 完
                break;
            }
            
            // ψ_1(r) 部分
            arb_sqr(s,Gamma_3,prec); // 1-（Gamma_3）^2 
            arb_neg(s,s);
            arb_add_si(s,s,1,prec);
            
            
            arb_mul(t,K_3_square,Gamma_3,prec);
            arb_neg(t,t);
            arb_add_si(t,t,1,prec);
            
            arb_div(t,t,s,prec); //系数完
            Help_psi_1_n(phi_1,r,order,prec);
            arb_mul(phi_1,t,phi_1,prec);
            
            // ∆ψ_1(r) 部分
            arb_div(t,K_3_square,Gamma_3,prec);
            arb_neg(t,t);
            arb_add_si(t,t,1,prec);
            arb_div(t,t,s,prec); //s用掉
            
            arb_sqr(s,R_3,prec);
            arb_div_si(s,s,3,prec);
            arb_mul(t,t,s,prec); //系数完
            
            Help_Laplacian_psi_1_n(D_phi_1,r,order,prec);
            arb_mul(D_phi_1,t,D_phi_1,prec);
            
            arb_add(s,phi_1,D_phi_1,prec); //右边完
            
            
            arb_mul(res,s,Mu_2,prec); // ζ(r) 完
            
            break;
        case power_law_type :
            
            //对应高斯型 ζ_G 的计算
            
            //包含两个参数 Mu_2、K_3_square, 这里，默认 ζ_∞ = 0
            //这里的求和分成两个部分：一是 ψ_1(r)，一是 ∆ψ_1(r)
            
            if( arb_equal(K_3_square,Gamma_3) )
            {
                Help_psi_1_n(phi_1,r,order,prec);
                
                arb_mul(res,phi_1,Mu_2,prec); // ζ(r) 完
                break;
            }
            printf("General -> typical_profile ->  zeta_Gauss_profile_n --> power_law_type 对于复杂的情况，还未实现");
            exit(1);
            break;
        case broken_power_law_type :
            
            //对应高斯型 ζ_G 的计算
            
            //包含两个参数 Mu_2、K_3_square, 这里，默认 ζ_∞ = 0
            //这里的求和分成两个部分：一是 ψ_1(r)，一是 ∆ψ_1(r)
            
            if( arb_equal(K_3_square,Gamma_3) )
            {
                Help_psi_1_n(phi_1,r,order,prec);
                
                arb_mul(res,phi_1,Mu_2,prec); // ζ(r) 完
                break;
            }
            printf("General -> typical_profile ->  zeta_Gauss_profile_n --> broken_power_law_type 对于复杂的情况，还未实现");
            exit(1);
            break;
        case box_type :
            //对应高斯型 ζ_G 的计算
            
            //包含两个参数 Mu_2、K_3_square, 这里，默认 ζ_∞ = 0
            //这里的求和分成两个部分：一是 ψ_1(r)，一是 ∆ψ_1(r)
            
            if( arb_equal(K_3_square,Gamma_3) )
            {
                Help_psi_1_n(phi_1,r,order,prec);
                
                arb_mul(res,phi_1,Mu_2,prec); // ζ(r) 完
                break;
            }
            printf("General -> typical_profile ->  zeta_Gauss_profile_n --> box_type 对于复杂的情况，还未实现");
            exit(1);
            break;
        case link_cmb_type :
            //对应高斯型 ζ_G 的计算
            
            //包含两个参数 Mu_2、K_3_square, 这里，默认 ζ_∞ = 0
            //这里的求和分成两个部分：一是 ψ_1(r)，一是 ∆ψ_1(r)
            
            if( arb_equal(K_3_square,Gamma_3) )
            {
                Help_psi_1_n(phi_1,r,order,prec);
                
                arb_mul(res,phi_1,Mu_2,prec); // ζ(r) 完
                break;
            }
            printf("General -> typical_profile ->  zeta_Gauss_profile_n --> link_cmb_type 对于复杂的情况，还未实现");
            exit(1);
            break;
        case delta_type :
            
            //对应高斯型 ζ_G 的计算
            
            //首先计算 x
            arb_mul(x,K_star,r,prec); //x=k_star*r
            
            switch(order)
            {
                case 0: //原函数
                    //功率谱为delta函数时，ζ(r) 易解析求出
                    // ζ(r)=mu_2*sinc(k_star*r)
                    
                    Help_sinc_n(s,x,0,prec);
                    
                    arb_mul(res,s,Mu_2,prec);
                    
                    break;
                    
                case 1: //一阶导
                    //功率谱为delta函数时易解析求出 x=k_star*r
                    //ζ(r)=mu_2*sinc(x)
                    //ζ(r)^prime=mu_2 * sinc'(x) * k_star
                    
                    Help_sinc_n(s,x,1,prec);
                    arb_mul(s,s,K_star,prec); //一阶导掉个 k_star
                    
                    arb_mul(res,s,Mu_2,prec);
                    
                    break;
                    
                case 2: //二阶导
                    //功率谱为delta函数时易解析求出 x=k_star*r
                    //ζ(r)=mu_2*sinc(x)
                    //ζ(r)^prime=mu_2 * 1/x *[cos(x)-sinc(x)] * k_star
                    //ζ(r)^prime^prime=mu_2 * (-1/x^2) * [x*sin(x)+2(cos(x)-sinc(x))] * (k_star)^2
                    
                    Help_sinc_n(s,x,2,prec);
                    arb_mul(s,s,K_star,prec); //二阶导掉个 (k_star)^2
                    arb_mul(s,s,K_star,prec);
                    
                    arb_mul(res,s,Mu_2,prec);
                    
                    break;
                    
                case 3: //三阶导
                    //ζ(r)=mu_2*sinc(x)  x=k_star*r
                    //求 ζ'''(r)
                    
                    Help_sinc_n(s,x,3,prec);
                    arb_mul(s,s,K_star,prec); //三阶导掉个 (k_star)^3
                    arb_mul(s,s,K_star,prec);
                    arb_mul(s,s,K_star,prec);
                    
                    arb_mul(res,s,Mu_2,prec);
                    
                    break;
                case 4: //四阶导
                    //ζ(r)=mu_2*sinc(x)   x=k_star*r
                    //求 ζ''''(r)
                    
                    Help_sinc_n(s,x,4,prec);
                    arb_mul(s,s,K_star,prec); //三阶导掉个 (k_star)^4
                    arb_mul(s,s,K_star,prec);
                    arb_mul(s,s,K_star,prec);
                    arb_mul(s,s,K_star,prec);
                    
                    arb_mul(res,s,Mu_2,prec);
                    
                    break;
                default:
                    printf("General -> typical_profile -> zeta_profile_n->delta_type->gaussian_type 输入有误\n");
                    exit(1);
            }
            
            
            break;
            
        default :
            printf("General -> typical_profile -> zeta_Gauss_profile_n 中的 Power_spectrum_type 有误\n");
            exit(1);
    }
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(x);
    
    arb_clear(phi_1);
    arb_clear(D_phi_1);
    
    return 0;
}


// ζ(r) 及其各阶导数
int zeta_profile_n(arb_t res, const arb_t r, const slong order, slong prec)
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
    //当 K_3_square=γ_3 时 ζ(r)=Mu_2 * ψ_1(r)
    // 而ζ(r)中的 ∆ψ(1) 部分，被 (R_3)^2 极大的压低，当 K_star 非常大时，基本可忽略
    
    //这里，不用对功率谱分类：lognormal_type和delta_type
    //直接统一，对 Zeta_type 分类即可
    
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
    
    return 0;
}


