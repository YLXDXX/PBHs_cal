#include "threshold.h"
#include <stdlib.h>
#include <arb_hypgeom.h> //特殊函数

//计算 q 参数，其为用来计算振幅的临界值
int Q_parameter(arb_t res, const arb_t r_m, slong prec)
{
    //函数中所用变量
    arb_t s,t,c;
    arb_init(s);
    arb_init(t);
    arb_init(c);
    
    //q=-[ C''*r^2 / (4*C*[1-C/f(w)]) ]
    //f(w)=3* (1+w)/(5+3w) 当w=1/3时，f(w)=2/3
    
    C_r_profile_n(s,r_m, 2, prec); //C''(r_m)
    
    //分子部分
    arb_sqr(t,r_m,prec); //C''*r^2
    arb_mul(s,s,t,prec);
    
    
    C_r_profile_n(c,r_m,0,prec); // C(r_m)
    
    //分母部分
    interior_assist_f_w(t,Equation_Of_State_w,prec); //f(w)
    
    arb_div(t,c,t,prec); //4*C*[1-C/f(w)]
    arb_neg(t,t);
    arb_add_si(t,t,1,prec);
    arb_mul(t,t,c,prec);
    arb_mul_si(t,t,4,prec);
    
    arb_div(s,s,t,prec);
    arb_neg(res,s); //最后取负号
    
    
    if(arb_is_negative(res))
    {
        printf("General -> threshold -> Q_parameter 错误：q_parameter为负 ");
        arb_printn(res, 15,0);printf("\n请检查各步计算\n");
        //arb_abs(i_q,q); //q不会取负值
        exit(1);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(c);
    
    return 0;     
}


//利用 q 参数，计算临界值 δ_c（q）
int Delta_c_q_parameter_simple(arb_t res, const arb_t q, slong prec)
{
    //函数中所用变量
    arb_t s,t,w,u,c_1,c_2,i_q;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(u);
    arb_init(c_1);
    arb_init(c_2);
    arb_init(i_q);
    
    arb_abs(i_q,q);
    
    arb_inv(c_1,i_q,prec); //1/q
    
    arb_mul_si(c_2,c_1,5,prec); //5/(2q)
    arb_div_si(c_2,c_2,2,prec);
    
    arb_neg(s,c_1); //前系数部分
    arb_exp(s,s,prec);
    arb_mul_si(s,s,4,prec);
    arb_div_si(s,s,15,prec);
    
    arb_neg(t,c_2); //后面分子部分
    arb_add_si(t,t,1,prec);
    arb_pow(t,i_q,t,prec);
    
    
    arb_gamma(w,c_2,prec);//Γ函数
    arb_hypgeom_gamma_upper(u,c_2,c_1,0,prec);//不完全Γ函数
    arb_sub(w,w,u,prec);
    
    
    
    arb_div(t,t,w,prec);
    arb_mul(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(u);
    arb_clear(c_1);
    arb_clear(c_2);
    arb_clear(i_q);
    
    return 0;  
}



static int interior_delta_c_q_parameter_new_F_1(arb_t res, const arb_t q, slong prec)
{
    arb_t a,b,c,d,s,t;
    
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(d);
    
    arb_init(s);
    arb_init(t);
    
    //利用了 Gauss hypergeometric function arb_hypgeom_2f1
    // 其4个参数，依次为a,b,c,d
    
    // 5/[2(1+q)]
    arb_add_ui(s,q,1,prec);
    arb_mul_ui(s,s,2,prec);
    arb_inv(s,s,prec);
    arb_mul_ui(s,s,5,prec);
    
    arb_one(a); //a=1
    
    arb_one(b); //b=1-5/[2(1+q)]
    arb_sub(b,b,s,prec);
    
    arb_one(c); //c=2-5/[2(1+q)]
    arb_mul_ui(c,c,2,prec);
    arb_sub(c,c,s,prec);
    
    arb_neg(d,q); //d=-q
    
    arb_hypgeom_2f1(t,a,b,c,d,0,prec);
    
    arb_set(res,t);
    
    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(d);
    arb_clear(s);
    arb_clear(t);
    
    return 0;
}

static int interior_delta_c_q_parameter_new_F_2(arb_t res, const arb_t q, const arb_t alpha, slong prec)
{
    arb_t a,b,c,d,s,t;
    
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(d);
    
    arb_init(s);
    arb_init(t);
    
    //利用了 Gauss hypergeometric function arb_hypgeom_2f1
    // 其4个参数，依次为a,b,c,d
    
    // 5/[2(1+q)]
    arb_add_ui(s,q,1,prec);
    arb_mul_ui(s,s,2,prec);
    arb_neg(t,s);  //后面用 t=-2(1+q)
    arb_inv(s,s,prec);
    arb_mul_ui(s,s,5,prec);
    
    arb_one(a); //a=1
    
    arb_one(b); //b=1-5/[2(1+q)]
    arb_sub(b,b,s,prec);
    
    arb_one(c); //c=2-5/[2(1+q)]
    arb_mul_ui(c,c,2,prec);
    arb_sub(c,c,s,prec);
    
    arb_one(d); // d = -q*(1-α)^[-2(1+q)]
    arb_sub(d,d,alpha,prec);
    arb_pow(d,d,t,prec); // t 用掉
    arb_mul(d,d,q,prec);
    arb_neg(d,d); //取负号
    
    
    arb_hypgeom_2f1(t,a,b,c,d,0,prec);
    
    arb_set(res,t);
    
    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(d);
    arb_clear(s);
    arb_clear(t);
    
    return 0;
}


static int interior_C_c_average_w(arb_t res, const arb_t w, slong prec)
{
    arb_t a,b,c,d,s,t;
    
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(d);
    
    arb_init(s);
    arb_init(t);
    
    //C_c(w)=a+b*arctan(c*w^d)
    
    arb_set_str(a,"-0.140381",prec);
    arb_set_str(b,"0.79538",prec);
    arb_set_str(c,"1.23593",prec);
    arb_set_str(d,"0.357491",prec);
    
    arb_pow(s,w,d,prec);
    arb_mul(s,s,c,prec);
    arb_atan(s,s,prec);
    
    arb_mul(s,s,b,prec);
    arb_add(s,s,a,prec);
    
    arb_set(res,s);
    
    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(d);
    arb_clear(s);
    arb_clear(t);
    return 0;
}



static int interior_q_parameter_alpha_w(arb_t res, const arb_t w, slong prec)
{
    arb_t e,f,g,h,s,t;
    
    arb_init(e);
    arb_init(f);
    arb_init(g);
    arb_init(h);
    
    arb_init(s);
    arb_init(t);
    
    //α(w)=e+f*arctan(g*w^h)
    
    arb_set_str(e,"2.00804",prec);
    arb_set_str(f,"-1.10936",prec);
    arb_set_str(g,"10.2801",prec);
    arb_set_str(h,"1.113",prec);
    
    
    arb_pow(s,w,h,prec);
    arb_mul(s,s,g,prec);
    arb_atan(s,s,prec);
    
    arb_mul(s,s,f,prec);
    arb_add(s,s,e,prec);
    
    arb_set(res,s);
    
    
    arb_clear(e);
    arb_clear(f);
    arb_clear(g);
    arb_clear(h);
    arb_clear(s);
    arb_clear(t);
    return 0;
}


static int interior_q_parameter_g_q_w(arb_t res, const arb_t q, const arb_t w, const arb_t alpha, slong prec)
{
    arb_t s,t;
    
    //arb_init(alpha);
    arb_init(s);
    arb_init(t);
    
    //g(q,w)= 3(1+q)/[ α(2q-3)*[3+α(α-3)] ]
    
    //求 α
    //interior_q_parameter_alpha_w(alpha,w,prec); 
    // α 直接作为参数传入，避免重复计算
    
    //分母
    arb_mul_ui(s,q,2,prec); //α(2q-3)
    arb_sub_ui(s,s,3,prec);
    arb_mul(s,s,alpha,prec);
    
    arb_sub_ui(t,alpha,3,prec); // [3+α(α-3)]
    arb_mul(t,t,alpha,prec);
    arb_add_ui(t,t,3,prec);
    
    arb_mul(s,s,t,prec);
    
    //分子
    arb_add_ui(t,q,1,prec); //3(1+q)
    arb_mul_ui(t,t,3,prec);
    
    arb_div(t,t,s,prec);
    
    arb_set(res,t);
    
    //arb_clear(alpha);
    arb_clear(s);
    arb_clear(t);
    return 0;
}

//新的解释方法，利用 q 参数，计算临界值 δ_c（q），参考 2007.05564
int Delta_c_q_parameter_new(arb_t res, const arb_t q, slong prec)
{
    
    arb_t alpha,w,s,t,u,z;
    
    arb_init(alpha);
    arb_init(w);
    arb_init(s);
    arb_init(t);
    arb_init(u);
    arb_init(z);
    
    //对于辐射，这里状态方程 w=1/3 
    //arb_one(w);
    //arb_div_ui(w,w,3,prec);
    arb_set(w,Equation_Of_State_w); //状态方程参数作为合局变量
    
    //求 α
    interior_q_parameter_alpha_w(alpha,w,prec);
    
    
    //前面系数 C_c(w)/g(q,w)
    
    interior_C_c_average_w(s, w, prec);//C_c(w)
    
    interior_q_parameter_g_q_w(t, q, w, alpha, prec);//g(q,w)
    
    arb_div(u,s,t,prec);
    
    
    //后面分母 -F_1(q)+(1-α)^(3-2q) * F_2(q,α)
    
    arb_one(s);//(1-α)
    arb_sub(s,s,alpha,prec);
    arb_mul_ui(t,q,2,prec);//(3-2q)
    arb_neg(t,t);
    arb_add_ui(t,t,3,prec);
    
    arb_pow(s,s,t,prec);
    
    interior_delta_c_q_parameter_new_F_2(t, q, alpha, prec); // F_2(q,α)
    interior_delta_c_q_parameter_new_F_1(z, q, prec); //F_1(q)
    
    arb_mul(s,s,t,prec);
    arb_neg(z,z);
    arb_add(z,z,s,prec);
    
    arb_inv(z,z,prec);
    arb_mul(u,u,z,prec);
    
    arb_set(res,u);
    
    arb_clear(alpha);
    arb_clear(w);
    arb_clear(s);
    arb_clear(t);
    arb_clear(u);
    arb_clear(z);
    return 0;
}


//面积半径R(r)
int Area_R(arb_t res, const arb_t r, slong prec)
{
    //函数中所用变量
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    //R(r)=a*r*exp[ζ(r)]
    
    arb_mul(t,Scale_factor_a,r,prec);//系数
    
    zeta_profile_n(s,r,0,prec);//指数
    arb_exp(s,s,prec);
    
    arb_mul(res,t,s,prec);
    
    arb_clear(s);
    arb_clear(t);
    
    return 0;
}


//求C_m的平均值用
//积分函数 C(r)* R(r)^2 * dR(r)
static int interior_C_m_average(arb_t res, const arb_t r, void* params, const slong order, slong prec)
{
    arb_t s,t,w;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //C(r)* R(r)^2 * dR(r) = C(r)*R(r)^3*(1/r +ζ')dr
    //R(r)=a*r*exp(zeta)
    //dR(r)=a*exp(ζ) * (1+r*ζ^preime)) * dr = R(r)* (1+r*ζ^preime))/r *dr =R(1/r +ζ')dr
    //其中 (1+r*ζ^preime))^2=1-3/2*C(r)
    
    
    //前半部分 C(r)*R^3(r)
    C_r_profile_n(t,r,0,prec);
    
    Area_R(s,r,prec);
    arb_pow_ui(s,s,3,prec); // R^3
    arb_mul(t,t,s,prec);
    
    
    
    //后部分 (1/r +ζ')dr
    arb_inv(s,r,prec);
    zeta_profile_n(w,r,1,prec); //一阶导
    arb_add(s,s,w,prec);
    
    arb_mul(res,t,s,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    
    return 0;
}


//求C_m的平均值
int C_m_average(arb_t res, const arb_t r_m, slong prec)
{
    //函数中所用变量
    arb_t s,a,b,alpha,V;
    
    arb_init(s);
    arb_init(a);
    arb_init(b);
    arb_init(alpha);
    arb_init(V);
    
    //积分区应[0,r_m]，2109.00791 误写为 [0, R_m]，「按照他的结果也可推得」
    //不过，若是[0, R_m]的话，就和尺度有关了
    //另外，r 的取值范围 [0,r_m] 刚合适
    
    //这里，积分的上下限采用不同的估算法式，结果不同
    //采用1907.13311的,[0, r_m]
    //而2007.05564作了更加细致的处理[r(α), r_m]
    
    //积分下限 a
    if(PT_Mu_th_METHOD==average_method_simple)
    {
        //a=0
        arb_set_str(a,"1E-100",prec); //设定值，不要设为零，里面有除以r的运算
        
    }else if(PT_Mu_th_METHOD==average_method_new)
    {
        
        
        arb_set(s,Equation_Of_State_w);
        interior_q_parameter_alpha_w(alpha,s,prec);
        
        //a=r(α)=r_m*(1-α)
        interior_q_parameter_alpha_w(alpha,s,prec);
        arb_neg(a,alpha);
        arb_add_ui(a,a,1,prec);
        arb_mul(a,a,r_m,prec);
        
    }
    
    
    //积分上限 b=r_m
    arb_set(b,r_m); 
    
    
    //∫_{r_min}^{r_max} C(r)R^2(r)dR(r)
    int ret_judge=0;
    ret_judge=Integration_arb(s, interior_C_m_average, NULL, 0,
                                        a, b, C_m_average_precision,
                                        C_m_average_iterate_min, C_m_average_iterate_max, prec );
    if(ret_judge==1)
    {
        printf("C_m_average  \t 达到最大迭代次数\n");
    }
    
    //求最后的平均值
    arb_mul_si(s,s,3,prec); //分子 3*int（...）
    
    Area_R(b,r_m,prec); // 分母 R^3
    arb_pow_ui(b,b,3,prec);
    
    arb_div(s,s,b,prec);
    
    //最后考虑到2007.05564的精细处理，还有个体积因子
    if(PT_Mu_th_METHOD==average_method_simple)
    {
        //V=1
        arb_set(res,s);
        
    }else if(PT_Mu_th_METHOD==average_method_new)
    {
        //V(α)=α[3+(α-3)α]
        
        arb_sub_ui(V,alpha,3,prec);
        arb_mul(V,V,alpha,prec);
        arb_add_ui(V,V,3,prec);
        arb_mul(V,V,alpha,prec);
        
        arb_div(res,s,V,prec);
        
    }
    
    
    
    arb_clear(s);
    arb_clear(a);
    arb_clear(b);
    arb_clear(alpha);
    arb_clear(V);
    
    return 0; 
}



//找PT_mu的临界值用 C(r)平均值版本
static int interior_PT_Mu_th_average(arb_t res, const arb_t mu, void * zeta_k, const slong order, slong prec)
{
    
    arb_t t,u,r_max,C_th_average;
    
    arb_init(t);
    arb_init(u);
    arb_init(r_max);
    arb_init(C_th_average);
    
    
    //需要修改 PT_mu，因为非高斯时 r_m 与 PT_mu 有关
    //当profile没有采取简化时，还与 PT_k 有关
    arb_set(PT_mu,mu);
    arb_set(PT_k,zeta_k);
    
    
    if(Stdout_verbose==true)
    {
        printf("\n当前尝试 μ_2 为： ");arb_printn(PT_mu, 15,0);printf("\n开始计算\nr_max： ");
    }
    
    Find_r_max(r_max, prec); //获得C(r)极大值的位置 r_max
    
    if(Stdout_verbose==true)
    {
        arb_printn(r_max, 15,0);printf("\nC(r_m):");
    }
    
    C_r_profile_n(t, r_max, 0,prec);
    
    if(Stdout_verbose==true)
    {
        arb_printn(t, 15,0);printf("\nC_m_average： ");
    }
    
    C_m_average(t,r_max,prec); //通过上面获得的r_m得到C_m_average
    
    if(Stdout_verbose==true)
    {
        arb_printn(t, 50,0);printf("\nC_th_average： ");
    }
    
    
    // mean compaction function threshold
    //采用1907.13311的 对于辐射来说w=1/3，C_th_average=2/5
    //而2007.05564作了更加细致的处理 C_th_average(w)
    if(PT_Mu_th_METHOD==average_method_simple)
    {
        //C_th_average=2/5
        arb_set_str(C_th_average,"0.4",prec);
        
    }else if(PT_Mu_th_METHOD==average_method_new)
    {
        //C_th_average=C_c(w)
        
        arb_set(u,Equation_Of_State_w);
        
        interior_C_c_average_w(C_th_average,u,prec);
    }
    
    if(Stdout_verbose==true)
    {
        arb_printn(C_th_average, 50,0);printf("\nC_m_average - C_th_average：");
    }
    
    arb_sub(res,t,C_th_average,prec);
    
    if(Stdout_verbose==true)
    {
        arb_printn(res, 15,0);printf("\n\n");
    }
    
    
    arb_clear(t);
    arb_clear(u);
    arb_clear(r_max);
    arb_clear(C_th_average);
    
    return 0; 
}


//找PT_mu的临界值用 q参数版本
static int interior_PT_Mu_th_q_parameter(arb_t res, const arb_t mu, void * zeta_k, const slong order, slong prec)
{
    
    arb_t s,t,r_max;
    
    arb_init(s);
    arb_init(t);
    arb_init(r_max);
    
    //需要修改 PT_mu，因为非高斯时 r_m 与 PT_mu 有关
    //当profile没有采取简化时，还与 PT_k 有关
    arb_set(PT_mu,mu);
    arb_set(PT_k,zeta_k);
    
    if(Stdout_verbose==true)
    {
        printf("\n当前尝试 μ_2 为： ");arb_printn(PT_mu, 15,0);printf("\n开始计算r_max: ");
    }
    
    
    Find_r_max(r_max, prec); //获得C(r)极大值的位置 r_max
    
    if(Stdout_verbose==true)
    {
        arb_printn(r_max, 15,0);printf("\nq_parameter: ");
    }
    
    
    Q_parameter(s,r_max,prec);; //通过上面获得的r_m得到 q_parameter
    
    if(Stdout_verbose==true)
    {
        arb_printn(s, 15,0);printf("\nδ_c(q): ");
    }
    
    
    //这里有两种方法
    if(PT_Mu_th_METHOD==q_parameter_method_simple)
    {
        Delta_c_q_parameter_simple(t,s,3*prec); // 由 q 得到 δ_c
        
    }else if(PT_Mu_th_METHOD==q_parameter_method_new)
    {
        Delta_c_q_parameter_new(t,s,3*prec); // 由 q 得到 δ_c
    }
    
    if(Stdout_verbose==true)
    {
        arb_printn(t, 50,0);printf("\nC(r_m)： ");
    }
    
    
    C_r_profile_n(s,r_max,0,prec); // C(r_max)
    
    if(Stdout_verbose==true)
    {
        arb_printn(s, 50,0);printf("\nδ_c(q) - C(r_m)： ");
    }
    
    
    arb_sub(res,t,s,prec);
    
    if(Stdout_verbose==true)
    {
        arb_printn(res, 15,0);printf("\n\n");
    }
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(r_max);
    
    return 0; 
}



//找PT_mu的临界值
int Find_PT_Mu_th(arb_t res, const arb_t zeta_k, slong prec)
{
    arb_t k;
    arb_init(k);
    
    arb_set(k,zeta_k);
    
    //求临界值有两种方法
    //求临界值的四种方式 q_parameter_method_simple / q_parameter_method_new
    //                average_method_simple / average_method_new
    if ( PT_Mu_th_METHOD==q_parameter_method_simple || PT_Mu_th_METHOD==q_parameter_method_new )
    {
        Find_interval_root(res, interior_PT_Mu_th_q_parameter, k, 0,
                           Int_mu_min, Int_mu_max, Int_mu_precision,
                           Root_mu_num, Root_Normal, prec);
    }else if ( PT_Mu_th_METHOD==average_method_simple || PT_Mu_th_METHOD==average_method_new )
    {
        Find_interval_root(res, interior_PT_Mu_th_average, k, 0,
                           Int_mu_min, Int_mu_max, Int_mu_precision,
                           Root_mu_num, Root_Normal, prec);
    }else
    {
        printf("General -> threshold -> Find_PT_Mu_th -> method 输入有误\n");
        exit(1);
    }
    
    printf("\n\nFind PT_mu_th: ");arb_printn(res, 50,0);printf("\n\n"); //打印变量
    
    arb_clear(k);
    
    return 0;
}


//将compaction func C 的值转为 C_l 的值
int Trans_C_to_C_l(arb_t res, const arb_t x, slong prec)
{
    arb_t s,t;
    
    arb_init(s);
    arb_init(t);
    
    //二次方程，有两个根，考虑到第一型扰动，取其中的一个根
    //C_l_th=4/3*( 1-sqrt(1-3/2*x) )
    
    arb_one(t); 
    arb_mul_ui(s,t,3,prec); //sqrt(1-3/2*x)
    arb_div_ui(s,s,2,prec);
    arb_mul(s,s,x,prec);
    arb_neg(s,s);
    arb_add_ui(s,s,1,prec);
    arb_sqrt(s,s,prec);
    
    arb_neg(s,s);
    arb_add_ui(s,s,1,prec);
    
    arb_mul_ui(s,s,4,prec);
    arb_div_ui(res,s,3,prec);
    
    arb_clear(s);
    arb_clear(t);
    
    return 0; 
}


static int interior_Get_PK_mu_max(arb_t res, const arb_t mu, void * zeta_k, const slong order, slong prec)
{
    arb_t s;
    arb_init(s);
    
    //与 μ 和 k 均有关
    arb_set(PT_mu,mu);
    arb_set(PT_k,zeta_k);
    
    //通过利用 C_l 的最大值 4/3 来得到 μ 的最大值
    //C_l=-4/3*r*ζ' --> r*ζ'的最大值为 -1
    
    zeta_profile_n(s, R_MAX, 1, prec);
    arb_mul(s,s,R_MAX,prec);
    
    arb_add_ui(res,s,1,prec);
    
    arb_clear(s);
    
    return 0; 
}

//求出参数 μ 的上限
int Get_PK_mu_max(arb_t res, const arb_t zeta_k, slong prec)
{
    arb_t s,a,b,k;
    
    arb_init(s);
    arb_init(a);
    arb_init(b);
    arb_init(k);
    
    //当不考虑profile的简化时，与k参数相关
    arb_set(k,zeta_k);
    
    //找根的最小值，从 PT_mu_th 开始
    arb_set(a,PT_mu_th); 
    arb_set(b,Int_mu_max);
    arb_mul_ui(b,b,2,prec);
    
    Find_interval_root(s, interior_Get_PK_mu_max, k, 0,
                       a, b, Int_mu_precision,
                       2*Root_mu_num, Root_Normal, prec);
    arb_set(res,s);
    
    arb_clear(s);
    arb_clear(a);
    arb_clear(b);
    arb_clear(k);
    
    return 0; 
}
