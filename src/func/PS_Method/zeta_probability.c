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
    //利用概率守恒求解 P(ζ)=P(ζ_G)*sqrt(1-h*zeta_G)
    //其中 ζ_G=1/h * (1- (1 - h/2 * ζ)^2 )
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
    arb_div_si(s,s,2,prec); //解析解需要注释掉这行
    arb_abs(s,s); //这里应该取绝对值，概率不应为负，这里本身=sqrt(1-h*R_G)
    
    
//     //归一化系数，解析求解 1+Erf(1/(h*sqrt(2*Σ_yy)))
//     arb_mul_ui(t,PS_Sigma_YY,2,prec);
//     arb_sqrt(t,t,prec);
//     arb_mul(t,t,h,prec);
//     arb_inv(t,t,prec);
//     arb_hypgeom_erf(t,t,prec);
//     arb_add_si(t,t,1,prec);
//     
//     arb_div(s,s,t,prec); //除以归一化系数
    
    
    arb_mul(res,w,s,prec);
    
    arb_clear(t);
    arb_clear(s);
    arb_clear(w);
    arb_clear(h);
    return 0;
}



static int interior_probability_up_step(arb_t res, const arb_t zeta_G, void* zeta, const slong order, slong prec)
{
    arb_t s;
    arb_init(s);
    
    Non_Gaussianity_up_step_n(s, zeta_G, 0, prec);
    arb_sub(res,s,zeta,prec);
    
    arb_clear(s);
    return 0;
}


static int interior_probability_narrow_1_up_step(arb_t res, const arb_t zeta_G, void* zeta, const slong order, slong prec)
{
    arb_t s;
    arb_init(s);
    
    Non_Gaussianity_narrow_1_up_step_n(s, zeta_G, 0, prec);
    arb_sub(res,s,zeta,prec);
    
    arb_clear(s);
    return 0;
}

static int interior_probability_narrow_1_2_up_step(arb_t res, const arb_t zeta_G, void* zeta, const slong order, slong prec)
{
    arb_t s;
    arb_init(s);
    
    Non_Gaussianity_narrow_1_2_up_step_n(s, zeta_G, 0, prec);
    arb_sub(res,s,zeta,prec);
    
    arb_clear(s);
    return 0;
}

//计算 ζ 的概率密度 P(ζ)
int Probability_zeta(arb_t res, const arb_t zeta, void* p, const slong order, slong prec)
{
    arb_t s,t,w,a,b,sum;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(a);
    arb_init(b);
    arb_init(sum);
    
    arb_zero(sum);
    
    int root_num; //根的个数
    arb_ptr muil_r; //存储多个根
    arb_ptr* m_r; //改变muil_r指针指向的地址，需要一个指向该指针的指针
    m_r=&muil_r;
    
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
            
            arb_set_ui(s,2);
            arb_div(s,s,Up_step_h,prec);
            arb_abs(s,s);
            
            if( arb_ge(zeta,s) ) //这里 ζ ≤ 2/|h|
            {
                arb_zero(res);
            }else
            {
            //利用概率守恒求解，P(ζ)= P(ζ_G)*dζ_G/dζ , 可能存在归一化的问题
            interior_probability_upward_step(s,zeta,NULL,0,prec); 
            arb_set(res,s); //解析求解归一化系数用
            }
            
            if(0)
            {
                //同样 是利用概率守恒求解，这里未利用解析表达式，纯数值计算
                //ζ=F(ζ_G)， 注意，这里的反函数 F^-1(ζ) 不是单值的，需要找出其所有的 ζ_G , 再求和
                //概率 P(ζ)= Σ P(ζ_G)*1/(dζ/dζ_G)
                //这里需要利用求根方法，找出所有的 ζ_G 根
                
                arb_set_ui(s,2);
                arb_div(s,s,Up_step_h,prec);
                arb_abs(s,s);
                
                if( arb_ge(zeta,s) ) //这里 ζ ≤ 2/|h|
                {
                    arb_zero(res);
                }else
                {
                    arb_set(w,zeta);
                    root_num=Find_interval_multi_root(m_r, interior_probability_up_step, w, 0,
                                                      PS_Root_zeta_to_zeta_G_min, PS_Root_zeta_to_zeta_G_max,
                                                      PS_Root_zeta_to_zeta_G_precision,
                                                      PS_Root_zeta_to_zeta_G_num,prec);
                    //printf("根的个数为: %i\n",root_num);
                    //将每个根对应的根率加起来
                    arb_zero(sum);
                    for(int root_i=0; root_i<root_num; root_i++)
                    {
                        //每个根对应的概率相加 P(ζ)= Σ P(ζ_G)*1/(dζ/dζ_G)
                        interior_probability_gauss(t, muil_r+root_i, prec);
                        Non_Gaussianity_up_step_n(s, muil_r+root_i, 1, prec);
                        arb_div(s,t,s,prec);
                        
                        arb_add(sum,sum,s,prec);
                    }
                    
                    arb_set(res,sum);
                    
                    _arb_vec_clear(muil_r, root_num); //清理数组
                }
            }
            
            break;
            
        case narrow_step_1_type :
            //ζ=F(ζ_G)， 注意，这里的反函数 F^-1(ζ) 不是单值的，需要找出其所有的 ζ_G , 再求和
            //概率 P(ζ)= Σ P(ζ_G)*1/(dζ/dζ_G)
            //这里需要利用求根方法，找出所有的 ζ_G 根
            
            arb_set(w,zeta);
            
            root_num=Find_interval_multi_root(m_r, interior_probability_narrow_1_up_step, w, 0,
                                              PS_Root_zeta_to_zeta_G_min, PS_Root_zeta_to_zeta_G_max,
                                              PS_Root_zeta_to_zeta_G_precision,
                                              PS_Root_zeta_to_zeta_G_num,prec);
            //printf("根的个数为: %i\n",root_num);
            //将每个根对应的根率加起来
            arb_zero(sum);
            for(int root_i=0; root_i<root_num; root_i++)
            {
                //每个根对应的概率相加 P(ζ)= Σ P(ζ_G)*1/(dζ/dζ_G)
                interior_probability_gauss(t, muil_r+root_i, prec);
                Non_Gaussianity_narrow_1_2_up_step_n(s, muil_r+root_i, 1, prec);
                arb_div(s,t,s,prec);
                
                arb_add(sum,sum,s,prec);
            }
            
            arb_set(res,sum);
            
            _arb_vec_clear(muil_r, root_num); //清理数组
            
        case narrow_step_1_2_type :
            //ζ=F(ζ_G)， 注意，这里的反函数 F^-1(ζ) 不是单值的，需要找出其所有的 ζ_G , 再求和
            //概率 P(ζ)= Σ P(ζ_G)*1/(dζ/dζ_G)
            //这里需要利用求根方法，找出所有的 ζ_G 根
            
            arb_set(w,zeta);
            
            root_num=Find_interval_multi_root(m_r, interior_probability_narrow_1_2_up_step, w, 0,
                                              PS_Root_zeta_to_zeta_G_min, PS_Root_zeta_to_zeta_G_max,
                                              PS_Root_zeta_to_zeta_G_precision,
                                              PS_Root_zeta_to_zeta_G_num,prec);
            //printf("根的个数为: %i\n",root_num);
            //将每个根对应的根率加起来
            arb_zero(sum);
            for(int root_i=0; root_i<root_num; root_i++)
            {
                //每个根对应的概率相加 P(ζ)= Σ P(ζ_G)*1/(dζ/dζ_G)
                interior_probability_gauss(t, muil_r+root_i, prec);
                Non_Gaussianity_narrow_1_2_up_step_n(s, muil_r+root_i, 1, prec);
                arb_div(s,t,s,prec);
                
                arb_add(sum,sum,s,prec);
            }
            
            arb_set(res,sum);
            
            _arb_vec_clear(muil_r, root_num); //清理数组
            
            break;
            
        default:
            printf("PS_Method -> zeta_probability -> Probability_zeta ->  Zeta_type 不正确\n" );
            exit(1);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(a);
    arb_clear(b);
    arb_clear(sum);
    
    return 0;
}


