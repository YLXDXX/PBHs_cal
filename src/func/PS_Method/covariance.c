#include "covariance.h"
#include <stdlib.h>

//计算协方差

//协方差的计算，跟功率谱有关，这里 delta_type 形式的功率谱
//采用如下形式 P(k)=A*δ[ln(k)-ln(k_star)]

//xx部分
int interior_variance_xx(arb_t res, const arb_t k, void* r, const slong order, slong prec)
{
    //局部变量设定
    arb_t t,s,w,x;
    
    //初始化参数
    arb_init(t);
    arb_init(s);
    arb_init(w);
    arb_init(x);
    
    //为计算方便，需作变量代换，取对数单位 ln(k)
    //其中 r 以参数形式传入
    // x^2 * [sinc'(x)]^2 * P(k)  其中 x= kr -> exp(k)*r
    //注意，这里的 sinc'(x) 跟前面不一样，前面是 d[sinc(x)]/dr，这里是 d[sinc(x)]/dx 直接求导，不会掉个k下来
    
    arb_exp(w,k,prec); //求 x=exp(k)*r //w后面用
    arb_mul(x,w,r,prec);
    
    //左边
    arb_sqr(s,x,prec);
    
    //中间 d[sinc(x)]/dx
    Help_sinc_n(t,x,1,prec); //一阶导
    arb_sqr(t,t,prec);
    arb_mul(s,s,t,prec);
    
    //右边功率谱部分
    power_spectrum(t, k, prec);
    
    
    //考虑转移函数
    if(Transfer_Function)
    {
        arb_mul(s,s,t,prec);
        
        //加入转移函数
        //arb_exp(w,k,prec); //注意，这里需作变量代换，取对数单位 ln(k)
        Linear_transfer_function(t,w,r,prec); //传入的w已恢复成原本的k
        arb_sqr(t,t,prec); //这里是转移函数的平方
        
        arb_mul(res,s,t,prec);
    }else
    {
        arb_mul(res,s,t,prec); //未加转移函数的结果
    }
    
    
    
    //完成计算，释放
    arb_clear(t);
    arb_clear(s);
    arb_clear(w);
    arb_clear(x);
    
    return 0;
}


//xx部分
int Variance_XX(arb_t res, const arb_t r, slong prec)
{
    //局部变量设定
    arb_t s,t,x,t_r;
    
    arb_init(s);
    arb_init(t);
    arb_init(x);
    arb_init(t_r);
    
    arb_set(t_r,r);
    
    switch(Power_spectrum_type)
    {
        case lognormal_type :
            
            //使用新的gauss_kronrod积分算法
            Integration_arb(res, interior_variance_xx, t_r, 0, 
                                      PS_Int_variance_min, PS_Int_variance_max,PS_Int_variance_precision,
                                      Integration_iterate_min,Integration_iterate_max, prec);
            
            break;
        case power_law_type :
            //使用新的gauss_kronrod积分算法
            Integration_arb(res, interior_variance_xx, t_r, 0, 
                                      PS_Int_variance_min, PS_Int_variance_max,PS_Int_variance_precision,
                                      Integration_iterate_min,Integration_iterate_max, prec);
            
            break;
        case broken_power_law_type :
            //使用新的gauss_kronrod积分算法
            Integration_arb(res, interior_variance_xx, t_r, 0, 
                                      PS_Int_variance_min, PS_Int_variance_max,PS_Int_variance_precision,
                                      Integration_iterate_min,Integration_iterate_max, prec);
            
            break;
        case box_type :
            //使用新的gauss_kronrod积分算法
            Integration_arb(res, interior_variance_xx, t_r, 0, 
                                      PS_Int_variance_min, PS_Int_variance_max,PS_Int_variance_precision,
                                      Integration_iterate_min,Integration_iterate_max, prec);
            break;
        case link_cmb_type :
            //使用新的gauss_kronrod积分算法
            Integration_arb(res, interior_variance_xx, t_r, 0, 
                                      PS_Int_variance_min, PS_Int_variance_max,PS_Int_variance_precision,
                                      Integration_iterate_min,Integration_iterate_max, prec);
            break;
        case delta_type :
            
            //为计算方便，需作变量代换，取对数单位 ln(k)
            //其中 r 以参数形式传入
            // Σ_xx = x^2 * [sinc'(x)]^2 * P(k)  其中 x= kr -> exp(k)*r
            //注意，这里的 sinc'(x) 跟前面不一样，前面是 d[sinc(x)]/dr，这里是 d[sinc(x)]/dx 直接求导，不会掉个k下来
            //在δ谱下 Σ_xx = (x)^2*[sinc'(x)]^2 * A 其中 x=K_star*r
            
            arb_mul(x,K_star,r,prec);
            
            //左边
            arb_sqr(s,x,prec);
            
            //右边 d[sinc(x)]/dx
            Help_sinc_n(t,x,1,prec); //一阶导
            arb_sqr(t,t,prec);
            arb_mul(s,s,t,prec);
            
            
            //考虑转移函数
            if(Transfer_Function)
            {
                arb_mul(s,s,Power_A,prec);
                
                //加入转移函数
                //arb_exp(w,k,prec); //注意，这里需作变量代换，取对数单位 ln(k)
                Linear_transfer_function(t,K_star,r,prec); //这里需要传入原本的未取对数的k
                arb_sqr(t,t,prec); //这里是转移函数的平方
                
                arb_mul(res,s,t,prec);
            }else
            {
                arb_mul(res,s,Power_A,prec); //未加转移函数的结果
            }
            
            
            break;
            
        default:
            printf("PS_Method -> covariance -> Variance_XX -> Power_spectrum_type->zeta_type 有误\n");
            exit(1);
    }
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(x);
    arb_clear(t_r);
    
    return 0;
}



//xy部分
int interior_variance_xy(arb_t res, const arb_t k, void* r, const slong order, slong prec)
{
    //局部变量设定
    arb_t t,s,w,x;
    
    //初始化参数
    arb_init(t);
    arb_init(s);
    arb_init(w);
    arb_init(x);
    
    //为计算方便，需作变量代换，取对数单位 ln(k)
    //其中 r 以参数形式传入
    // x * sinc(x) * sinc'(x) * P(k)  其中 x= kr -> exp(k)*r
    //注意，这里的 sinc'(x) 跟前面不一样，前面是 d[sinc(x)]/dr，这里是 d[sinc(x)]/dx 直接求导，不会掉个k下来
    
    arb_exp(w,k,prec); //求 x=exp(k)*r //w后面用
    arb_mul(x,w,r,prec);
    
    //左边 x * sinc(x) * sinc'(x)
    Help_sinc_n(s,x,0,prec); //没求导
    Help_sinc_n(t,x,1,prec); //一阶导
    
    arb_mul(s,s,t,prec);
    arb_mul(s,s,x,prec); 
    
    
    //右边功率谱部分
    power_spectrum(t, k, prec);
    
    
    //考虑转移函数
    if(Transfer_Function)
    {
        arb_mul(s,s,t,prec);
        
        //加入转移函数
        //arb_exp(w,k,prec); //注意，这里需作变量代换，取对数单位 ln(k)
        Linear_transfer_function(t,w,r,prec); //传入的w已恢复成原本的k
        arb_sqr(t,t,prec); //这里是转移函数的平方
        
        arb_mul(res,s,t,prec);
    }else
    {
        arb_mul(res,s,t,prec); //未加转移函数的结果
    }
    
    //完成计算，释放
    arb_clear(t);
    arb_clear(s);
    arb_clear(w);
    arb_clear(x);
    
    return 0;
}
//xy部分
int Variance_XY(arb_t res, const arb_t r, slong prec)
{
    //局部变量设定
    arb_t s,t,x,t_r;
    
    arb_init(s);
    arb_init(t);
    arb_init(x);
    arb_init(t_r);
    
    arb_set(t_r,r);
    
    
    switch(Power_spectrum_type)
    {
        case lognormal_type :
            //使用新的gauss_kronrod积分算法
            Integration_arb(res, interior_variance_xy, t_r, 0, 
                                      PS_Int_variance_min, PS_Int_variance_max,PS_Int_variance_precision,
                                      Integration_iterate_min,Integration_iterate_max, prec);
            
            break;
        case power_law_type :
            //使用新的gauss_kronrod积分算法
            Integration_arb(res, interior_variance_xy, t_r, 0, 
                                      PS_Int_variance_min, PS_Int_variance_max,PS_Int_variance_precision,
                                      Integration_iterate_min,Integration_iterate_max, prec);
            
            break;
        case broken_power_law_type :
            //使用新的gauss_kronrod积分算法
            Integration_arb(res, interior_variance_xy, t_r, 0, 
                                      PS_Int_variance_min, PS_Int_variance_max,PS_Int_variance_precision,
                                      Integration_iterate_min,Integration_iterate_max, prec);
            break;
        case box_type :
            //使用新的gauss_kronrod积分算法
            Integration_arb(res, interior_variance_xy, t_r, 0, 
                                      PS_Int_variance_min, PS_Int_variance_max,PS_Int_variance_precision,
                                      Integration_iterate_min,Integration_iterate_max, prec);
            
            break;
        case link_cmb_type :
            //使用新的gauss_kronrod积分算法
            Integration_arb(res, interior_variance_xy, t_r, 0, 
                                      PS_Int_variance_min, PS_Int_variance_max,PS_Int_variance_precision,
                                      Integration_iterate_min,Integration_iterate_max, prec);
            break;
        case delta_type :
            
            //为计算方便，需作变量代换，取对数单位 ln(k)
            //其中 r 以参数形式传入
            // Σ_xy = x * sinc(x) * sinc'(x) * P(k)  其中 x= kr -> exp(k)*r
            //注意，这里的 sinc'(x) 跟前面不一样，前面是 d[sinc(x)]/dr，这里是 d[sinc(x)]/dx 直接求导，不会掉个k下来
            //在δ谱下 Σ_xy = x * sinc(x) * sinc'(x) * A 其中 x=K_star*r
            
            arb_mul(x,K_star,r,prec);
            
            // x * sinc(x) * sinc'(x)
            Help_sinc_n(s,x,0,prec); //没求导
            Help_sinc_n(t,x,1,prec); //一阶导
            
            arb_mul(s,s,t,prec);
            arb_mul(s,s,x,prec); 
            
            
            //考虑转移函数
            if(Transfer_Function)
            {
                arb_mul(s,s,Power_A,prec);;
                
                //加入转移函数
                //arb_exp(w,k,prec); //注意，这里需作变量代换，取对数单位 ln(k)
                Linear_transfer_function(t,K_star,r,prec); //这里需要传入原本的未取对数的k
                arb_sqr(t,t,prec); //这里是转移函数的平方
                
                arb_mul(res,s,t,prec);
            }else
            {
                arb_mul(res,s,Power_A,prec); //未加转移函数的结果
            }
            
            break;
            
        default:
            printf("PS_Method -> covariance -> Variance_XY -> Power_spectrum_type->zeta_type 有误\n");
            exit(1);
    }
    
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(x);
    arb_clear(t_r);
    
    return 0;
}

//yy部分
int interior_variance_yy(arb_t res, const arb_t k, void* r, const slong order, slong prec)
{
    //局部变量设定
    arb_t t,s,w,x;
    
    //初始化参数
    arb_init(t);
    arb_init(s);
    arb_init(w);
    arb_init(x);
    
    //为计算方便，需作变量代换，取对数单位 ln(k)
    //其中 r 以参数形式传入
    // [sinc(x)]^2 * P(k)  其中 x= kr -> exp(k)*r
    
    arb_exp(w,k,prec); //求 x=exp(k)*r //w后面用
    arb_mul(x,w,r,prec);
    
    //左边 [sinc(x)]^2
    Help_sinc_n(s,x,0,prec); //没求导
    arb_sqr(s,s,prec);
    
    //右边功率谱部分
    power_spectrum(t, k, prec);
    
    
    //考虑转移函数
    if(Transfer_Function)
    {
        arb_mul(s,s,t,prec);
        
        //加入转移函数
        //arb_exp(w,k,prec); //注意，这里需作变量代换，取对数单位 ln(k)
        Linear_transfer_function(t,w,r,prec); //传入的w已恢复成原本的k
        arb_sqr(t,t,prec); //这里是转移函数的平方
        
        arb_mul(res,s,t,prec);
    }else
    {
        arb_mul(res,s,t,prec); //未加转移函数的结果
    }
    
    
    //完成计算，释放
    arb_clear(t);
    arb_clear(s);
    arb_clear(w);
    arb_clear(x);
    
    return 0;
}
//yy部分
int Variance_YY(arb_t res, const arb_t r, slong prec)
{
    //局部变量设定
    arb_t s,t,x,t_r;
    
    arb_init(s);
    arb_init(t);
    arb_init(x);
    arb_init(t_r);
    
    arb_set(t_r,r);
    
    switch(Power_spectrum_type)
    {
        case lognormal_type :
            
            //使用新的gauss_kronrod积分算法
            Integration_arb(res, interior_variance_yy, t_r, 0, 
                                      PS_Int_variance_min, PS_Int_variance_max,PS_Int_variance_precision,
                                      Integration_iterate_min,Integration_iterate_max, prec);
            
            break;
        case power_law_type :
            //使用新的gauss_kronrod积分算法
            Integration_arb(res, interior_variance_yy, t_r, 0, 
                                      PS_Int_variance_min, PS_Int_variance_max,PS_Int_variance_precision,
                                      Integration_iterate_min,Integration_iterate_max, prec);
            
            break;
        case broken_power_law_type :
            //使用新的gauss_kronrod积分算法
            Integration_arb(res, interior_variance_yy, t_r, 0, 
                                      PS_Int_variance_min, PS_Int_variance_max,PS_Int_variance_precision,
                                      Integration_iterate_min,Integration_iterate_max, prec);
            
            break;
        case box_type :
            //使用新的gauss_kronrod积分算法
            Integration_arb(res, interior_variance_yy, t_r, 0, 
                                      PS_Int_variance_min, PS_Int_variance_max,PS_Int_variance_precision,
                                      Integration_iterate_min,Integration_iterate_max, prec);
            
            break;
        case link_cmb_type :
            //使用新的gauss_kronrod积分算法
            Integration_arb(res, interior_variance_yy, t_r, 0, 
                                      PS_Int_variance_min, PS_Int_variance_max,PS_Int_variance_precision,
                                      Integration_iterate_min,Integration_iterate_max, prec);
            
            break;
        case delta_type :
            
            //为计算方便，需作变量代换，取对数单位 ln(k)
            //其中 r 以参数形式传入
            //Σ_yy = [sinc(x)]^2 * P(k)  其中 x= kr -> exp(k)*r
            //在δ谱下 Σ_yy =  [sinc(x)]^2 * A 其中 x= K_star*r
            
            arb_mul(x,K_star,r,prec);
            
            // [sinc(x)]^2
            Help_sinc_n(s,x,0,prec); //没求导
            arb_sqr(s,s,prec);
            
            
            //考虑转移函数
            if(Transfer_Function)
            {
                arb_mul(s,s,Power_A,prec);
                
                //加入转移函数
                //arb_exp(w,k,prec); //注意，这里需作变量代换，取对数单位 ln(k)
                Linear_transfer_function(t,K_star,r,prec); //这里需要传入原本的未取对数的k
                arb_sqr(t,t,prec); //这里是转移函数的平方
                
                arb_mul(res,s,t,prec);
            }else
            {
                arb_mul(res,s,Power_A,prec); //未加转移函数的结果
            }
            
            break;
            
        default:
            printf("PS_Method -> covariance -> Variance_YY -> Power_spectrum_type->zeta_type 有误\n");
            exit(1);
    }
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(x);
    arb_clear(t_r);
    
    return 0;
}




 
