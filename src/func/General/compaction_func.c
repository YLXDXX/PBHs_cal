#include "compaction_func.h"
#include <stdlib.h>

void interior_assist_f_w(arb_t res, const arb_t w, slong prec)
{
    //f(w)=3*(1+w)/(5+3w)
    arb_t s,t;
    
    arb_init(s);
    arb_init(t);
    
    arb_add_ui(s,w,1,prec);//3*(1+w)
    arb_mul_ui(s,s,3,prec);
    
    arb_mul_ui(t,w,3,prec);
    arb_add_ui(t,t,5,prec);
    
    arb_div(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
}


// C(r) 及其各阶导数
int C_r_profile_n(arb_t res, const arb_t r, const slong order, slong prec)
{
    //函数中所用变量
    arb_t s,t,w,z_1,z_2,z_3;
    
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(z_1);
    arb_init(z_2);
    arb_init(z_3);
    
    //据曲率 ζ 类型，来得出相应的C(r)
    //C(r)=f(w) * [ 1 - (1 + r*ζ(r)^prime)^2 ]    其中 f(w)=3* (1+w)/(5+3w)
    //当w=1/3时，f(w)=2/3
    
    switch(order)
    {
        case 0: //原函数
            //C(r)=-f(w) * [ 2r*ζ' + (r*ζ')^2 ]
            zeta_profile_n(z_1,r,1,prec);
            
            arb_mul(s,z_1,r,prec);
            arb_mul_si(t,s,2,prec);
            
            arb_sqr(s,s,prec);
            arb_add(t,t,s,prec);
            
            interior_assist_f_w(s,Equation_Of_State_w,prec); //f(w)
            arb_mul(t,t,s,prec);
            
            arb_neg(res,t); //最后取负号
            
            break;
            
        case 1: //一阶导数 
            //C'(r)= -f(w)*2*(1+r*ζ')*(ζ'+r*ζ'')
            
            zeta_profile_n(z_1,r,1,prec);
            zeta_profile_n(z_2,r,2,prec);
            
            interior_assist_f_w(t,Equation_Of_State_w,prec); //f(w)*2
            arb_mul_ui(t,t,2,prec);
            
            arb_mul(s,r,z_1,prec); //-f(w)*2*(1+r*ζ')
            arb_add_si(s,s,1,prec);
            arb_mul(s,s,t,prec);
            arb_neg(s,s); //取负号
            
            arb_mul(t,r,z_2,prec);
            arb_add(t,t,z_1,prec);
            
            arb_mul(res,s,t,prec);
            
            break;
            
        case 2: //二阶导数
            //C''(r)/(-f(w))=4*ζ''+2*r*ζ'''+2*[ζ']^2+8*r*ζ'*ζ''+2*r^2*[ζ'']^2+2*r^2*ζ'*ζ'''
            //化简 C''(r) 可得
            // C''(r)=-f(w)*2*[ (r*ζ'+1)*(r*ζ'''+2*ζ'') + (r*ζ''+ζ')^2 ]
            
            zeta_profile_n(z_1,r,1,prec);
            zeta_profile_n(z_2,r,2,prec);
            zeta_profile_n(z_3,r,3,prec);
            
            //括号中前半部分
            arb_mul(s,r,z_1,prec);
            arb_add_si(s,s,1,prec);
            
            arb_mul(t,r,z_3,prec);
            arb_mul_si(w,z_2,2,prec);
            arb_add(t,t,w,prec);
            arb_mul(s,s,t,prec);
            
            //括号中后半部分
            arb_mul(w,r,z_2,prec);
            arb_add(t,w,z_1,prec);
            arb_sqr(t,t,prec);
            
            arb_add(s,s,t,prec);//括号完
            
            
            interior_assist_f_w(t,Equation_Of_State_w,prec); //f(w)*2
            arb_mul_ui(t,t,2,prec);
            
            arb_mul(s,s,t,prec);
            arb_neg(res,s); //最后取负号
            
            
            break;
            
        default:
            printf("General -> compaction_func -> C_r_profile_n 中 order 输入有误\n");
            exit(1);
    }
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(z_1);
    arb_clear(z_2);
    arb_clear(z_3);
    
    return 0;  
}



//求C(r)的最大值用
int interior_find_r_max_general(arb_t res, const arb_t r,
                                       void * param, const slong order,  slong prec)
{
    arb_t t,s;
    arb_init(t);
    arb_init(s);
    
    //d[C(r)]/dr=0 --> ζ’ + rζ’’=0
    
    //C_r_profile_n(res,r,1,prec); //采用最初的形式
    //return 0;
    
    //采用 ζ’ + rζ’’ 模式
    zeta_profile_n(t,r,1,prec);
    
    zeta_profile_n(s,r,2,prec);
    
    arb_mul(s,s,r,prec);
    
    arb_add(res,t,s,prec);
    
    
    arb_clear(t);
    arb_clear(s);
    
    return 0;
}



//求C(r)的最大值，在区间[a, b]中找
int Find_r_max(arb_t res, slong prec)
{
    //注意，这里求最大值要用导数来求零点，将Find_interval_root参数设数1
    //这样才能找出最大值的点，对于有些特殊的profile，第一个极值点可能是最小值点
    
    
    arb_t t;
    arb_init(t);
    
    //都采用 ζ’ + rζ’’=0 求根
    Find_interval_root(res, interior_find_r_max_general,NULL,0,
                       Int_r_min, Int_r_max, Int_r_precision,
                       Root_r_num, Root_C_Max, prec);
    arb_clear(t);
    return 0;
    
    //功率谱类型判断
    switch(Power_spectrum_type) 
    {
        case lognormal_type :
            
            //d[C(r)]/dr=0 --> ζ’ + rζ’’=0
            
            Find_interval_root(res, interior_find_r_max_general,NULL,0,
                               Int_r_min, Int_r_max, Int_r_precision,
                               Root_r_num, Root_C_Max, prec);
            
            break;
        case delta_type:
            //曲率 ζ 类型判断
            switch(Zeta_type) 
            {
                //注意，这里是通过找根的方式来找r_m
                //所给函数的概可能有很多个，需要找的是，正的且离零点最近的一个根???
                
                case gaussian_type :
                    //功率谱为delta函数时易解析求出
                    //C(r)取最值时满足的方程如下 x=k_star*r
                    // -k_star * (1/x^2) * ( x*cos(x) + (x^2-1)*sin(x) ) =0
                    //可以看到极值与PT_mu无关
                    
                    //基于解析的计算，最终可化为如下函数求解
                    // x*cos(x)+(pow(x,2)-1)*sin(x)=0
                    
                    //利用arb直接用求根法求解
                    Find_interval_root(res, interior_find_r_max_general, NULL, 0,
                                       Int_r_min, Int_r_max, Int_r_precision,
                                       Root_r_num, Root_C_Max, prec);
                    
                    break;
                    
                case power_expansion_type :
                    //功率谱为delta函数时易解析求出
                    //C(r)取最值时满足的方程如下 x=k_star*r
                    //5*x*[ sin(x)-x(cos(x)+x*sinx) ] - mu_2 * f_NL*[ 6+6*(x^2 - 1)*cos(2*x)-9*sin(2*x) ]
                    
                    //注意，这里的极值与mu_2有关
                    
                    //功率谱为delta函数时易解析求出
                    //C(r)取最值时满足的方程如下 x=k_star*r
                    //论文中的公式有问题
                    //5*x*[ sin(x)-x(cos(x)+x*sinx) ] - mu_2 * f_NL*[ 6+6*(x^2 - 1)*cos(2*x)-9*sin(2*x) ]
                    //5*x*[ sin(x)-x*cos(x)-x^2*sin(x) ] - mu_2 * f_NL*[ 6+6*(x^2 - 1)*cos(2*x)-9*sin(2*x) ]
                    //正确公式应为：
                    //5*x*[ sin(x)-x*cos(x)-x^2*sin(x)]+3*f_Nl*mu_2*[ x^2*(1+cos(2x)-sin(2x))-3*x*sin(2x)-2*cos(2x)+2 ]
                    
                    //利用arb直接用求根法求解
                    Find_interval_root(res, interior_find_r_max_general, NULL, 0,
                                       Int_r_min, Int_r_max, Int_r_precision,
                                       Root_r_num, Root_C_Max, prec);
                    
                    //r_m_func_f_NL(t,t,prec);
                    //arb_printn(t, 50,0);printf("\n"); //打印变量
                    //arb_div(res,res,K_star,prec); //前面的结果是 x=k_star*r
                    
                    break;
                default :
                    printf("General -> compaction_func -> Find_r_max -> zeta_type 有误\n");
                    exit(1);
            }
            
            break;
        default:
            printf("General -> compaction_func -> Find_r_max-> Power_spectrum_type 输入有误\n");
            exit(1);
    }
    
    return 0;     
}


 
