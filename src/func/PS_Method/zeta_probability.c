#include "zeta_probability.h"
#include <stdlib.h>
#include <arb_hypgeom.h> //误差函数用

//各种 ζ 类型的概率密度

//高斯概率密度 P_G=1/sqrt(2*Pi*Σ_YY)*exp[-(ζ_G)^2/(2*Σ_YY)]
int interior_probability_gauss(arb_t res, const arb_t zeta, slong prec)
{
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    //指数部分
    arb_sqr(s,zeta,prec);
    arb_mul_si(t,PS_Sigma_YY,2,prec);
    arb_div(s,s,t,prec);
    arb_neg(s,s);
    arb_exp(s,s,prec);
    
    //系数部分
    arb_mul(t,PS_Sigma_YY,Pi_2,prec);
    arb_sqrt(t,t,prec);
    arb_inv(t,t,prec);
    
    arb_mul(res,s,t,prec);
    
    //完成计算，释放
    arb_clear(s);
    arb_clear(t);
    
    return 0;
}



//exponential_tail 概率密度 ζ=-1/3 * ln(1-3*ζ_G)
int interior_probability_exponential_tail(arb_t res, const arb_t zeta, slong prec)
{  
    arb_t s,t,w;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //指数部分
    arb_mul_ui(w,zeta,3,prec);
    arb_neg(w,w);
    arb_set(t,w); //t后面指数上减法用
    arb_exp(w,w,prec);
    arb_sub_si(w,w,1,prec); //因后面要平方，相减顺序不重要
    arb_sqr(w,w,prec);
    
    arb_mul_si(s,PS_Sigma_YY,18,prec);
    arb_div(w,w,s,prec);
    arb_neg(w,w);
    arb_add(w,w,t,prec);
    
    arb_exp(w,w,prec);
    
    //系数部分
    arb_mul(t,PS_Sigma_YY,Pi_2,prec);
    arb_sqrt(t,t,prec);
    arb_inv(t,t,prec); //Sets z to 1/𝑥
    
    arb_mul(res,t,w,prec);
    
    //完成计算，释放
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    
    return 0;
}


//upward_step的概率分布  R=-2/|h|*[sqrt(1-|h|R_G)-1]
int interior_probability_upward_step(arb_t res, const arb_t zeta, void* p, const slong order, slong prec)
{
    //这里未归一化，可以解析求解归一化系数，亦可数值积分求归一化系数
    //这里解析求解积分有点难算，用数值积分
    //其中P(ζ)=P(ζ_G)*sqrt(1-h*zeta_G)，而ζ_G=1/h * (1- (1 - h/2 * ζ)^2 )
    arb_t t,s,w,h;
    
    arb_init(t);
    arb_init(s);
    arb_init(w);
    arb_init(h);
    
    //先通过ζ求出ζ_G 
    //可通过表达示化简求
    //R_G=R/4 * (4-h*R)
    
    arb_abs(h,Up_step_h);
    
    arb_mul(t,h,zeta,prec);
    arb_neg(t,t);
    arb_add_si(t,t,4,prec);
    arb_mul(t,t,zeta,prec);
    arb_div_si(t,t,4,prec); //t=R_G
    
    //再通过 ζ_G 求 P(ζ)
    //可通过表达示化简
    //Gauss(R_G)*sqrt(1-h*R_G)=Gauss(R_G) * (2-h*R)/2
    
    interior_probability_gauss(w,t,prec); //利用高斯的概率密度
    
    arb_mul(s,h,zeta,prec); //(2-h*R)/2
    arb_neg(s,s);
    arb_add_si(s,s,2,prec);
    //arb_div_si(s,s,2,prec); //解析解需要注释掉这行
    arb_abs(s,s); //这里应该取绝对值，概率不应为负，这里本身=sqrt(1-h*R_G)
    
    
    //归一化系数，解析求解 1+Erf(1/(h*sqrt(2*Σ_yy)))
    arb_mul_ui(t,PS_Sigma_YY,2,prec);
    arb_sqrt(t,t,prec);
    arb_mul(t,t,h,prec);
    arb_inv(t,t,prec);
    arb_hypgeom_erf(t,t,prec);
    arb_add_si(t,t,1,prec);
    
    
    arb_div(s,s,t,prec); //除以归一化系数
    
    
    arb_mul(res,w,s,prec);
    
    arb_clear(t);
    arb_clear(s);
    arb_clear(w);
    arb_clear(h);
    return 0;
}


//计算 ζ 的概率密度 P(ζ)
int Probability_zeta(arb_t res, const arb_t zeta, void* p, const slong order, slong prec)
{
    arb_t s,a,b;
    arb_init(s);
    arb_init(a);
    arb_init(b);
    
    //针对各种不同的 ζ 类型 分开讨论
    switch(Zeta_type) 
    {
        case gaussian_type :
            
            interior_probability_gauss(s,zeta,prec);
            arb_set(res,s);
            break;
            
        case exponential_tail_type :
            
            interior_probability_exponential_tail(s,zeta,prec);
            arb_set(res,s);
            break;
            
        case up_step_type :
            
            
            //求归一化系数
            //需要归一化，用数值积分来归一化
            //积分区间为 [-infinity，1/h] -> [a,b]
            
            //注意，ζ<2/h，ζ_G<1/h 两者的取值范围不一样
            
            /*
            arb_set_str(a,"-1.5",prec); //这里 -infinity 取1
            
            arb_abs(b,Up_step_h);
            arb_inv(b,b,prec);
            arb_mul_si(b,b,2,prec);
            
            arb_set_str(s,"1E-30",prec); //精度
            
            if(P_normalization_coefficient==NULL)
            {
                P_normalization_coefficient=_arb_vec_init(1); //只需要1个点
                //使用新的gauss_kronrod积分算法
                integration_gauss_kronrod(P_normalization_coefficient, interior_probability_upward_step, NULL, 0, 
                                          a, b,s,
                                          Integration_iterate_min,Integration_iterate_max, prec);
            }
            */
            
            interior_probability_upward_step(s,zeta,NULL,0,prec);
            
            //arb_div(res,s,P_normalization_coefficient,prec); //积分求解归一化系数用
            
            arb_set(res,s); //解析求解归一化系数用
            
            break;
            
        default:
            printf("PS_Method -> zeta_probability -> Probability_zeta ->  Zeta_type 不正确\n" );
            exit(1);
    }
    
    arb_clear(s);
    arb_clear(a);
    arb_clear(b);
    return 0;
}


