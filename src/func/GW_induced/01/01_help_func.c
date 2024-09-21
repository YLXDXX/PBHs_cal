#include "01_power_spectra.h"
#include <stdlib.h>
#include <arb_hypgeom.h>

//辐射时代的转移函数
static void Radiation_transfer_function(arb_t res, const arb_t z, void* param, const slong order, slong prec)
{
    arb_t s,t,w,sqrt_3;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(sqrt_3);
    
    arb_one(s);// s=z/sqrt{3}
    arb_mul_ui(s,s,3,prec);
    arb_sqrt(sqrt_3,s,prec);
    arb_div(s,z,sqrt_3,prec);
    
    switch(order) 
    {
        case 0: //原函数
            // 9/z^2 * ( sin(s)/s - cos(s) ) = 3/s^2 * ( sin(s)/s - cos(s) )
            arb_sin(t,s,prec); //sin(s)/s - cos(s)
            arb_div(t,t,s,prec);
            arb_cos(w,s,prec);
            arb_sub(t,t,w,prec);
            
            arb_mul_ui(t,t,3,prec);
            arb_sqr(w,s,prec);
            arb_div(res,t,w,prec);
            
            break;
        case 1: //一阶导
            // 1/sqrt{3} * ( (3*s^2-9)*sin(s)+9*s*cos(s) )/s^4
            arb_sqr(t,s,prec); //(3*s^2-9)*sin(s)
            arb_mul_ui(t,t,3,prec);
            arb_sub_ui(t,t,9,prec);
            arb_sin(w,s,prec);
            arb_mul(t,t,w,prec);
            
            arb_cos(w,s,prec); //9*s*cos(s)
            arb_mul(w,w,s,prec);
            arb_mul_ui(w,w,9,prec);
            
            arb_add(t,t,w,prec);
            arb_pow_ui(w,s,4,prec);
            arb_div(t,t,w,prec);
            
            arb_div(res,t,sqrt_3,prec);
            
            break;
        default:
            printf(" radiation_transfer_function 输入有误\n");
            exit(1);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(sqrt_3);
}

//
//数值积分部分
//

//此处算诱导引力波的方法来自：1804.07732
static int interior_GW_I_c_func(arb_t res, const arb_t tau, void* param, const slong order, slong prec)
{
    arb_t s,t,w,x_tau,y_tau,T_x_tau,T_y_tau,T_x_tau_prime,T_y_tau_prime;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(x_tau);
    arb_init(y_tau);
    arb_init(T_x_tau);
    arb_init(T_y_tau);
    arb_init(T_x_tau_prime);
    arb_init(T_y_tau_prime);
    
    //接收传入参数
    struct Func_transfer_parameter *func_x_y; //这里不需要手动分配，只需将其指向传入的指针即可
    func_x_y=param;
    
    arb_mul(x_tau,func_x_y->p_1,tau, prec); //x*τ
    arb_mul(y_tau,func_x_y->p_2,tau, prec); //y*τ
    
    Radiation_transfer_function(T_x_tau,x_tau,NULL,0,prec); //T(x*τ)
    Radiation_transfer_function(T_y_tau,y_tau,NULL,0,prec); //T(y*τ)
    
    Radiation_transfer_function(T_x_tau_prime,x_tau,NULL,1,prec); //T'(x*τ)
    Radiation_transfer_function(T_y_tau_prime,y_tau,NULL,1,prec); //T'(y*τ)
    
    
    //积分函数形式为$ \tau( - \sin\tau)\cdot4\Big\{2\mathcal{T}(x\tau)\mathcal{T}(y\tau) + \Big[\mathcal{T}(x\tau) + x\tau \mathcal{T}'(x\tau)\Big]\Big[\mathcal{T}(y\tau) + y\tau \mathcal{T}'(y\tau)\Big]\Big\} $
    //右边 I_c 与 I_s 共同部分
    arb_mul(s,x_tau,T_x_tau_prime,prec);
    arb_add(s,s,T_x_tau,prec);
    
    arb_mul(t,y_tau,T_y_tau_prime,prec);
    arb_add(t,t,T_y_tau,prec);
    arb_mul(s,s,t,prec);
    
    arb_mul(t,T_x_tau,T_y_tau,prec);
    arb_mul_ui(t,t,2,prec);
    arb_add(s,s,t,prec);
    arb_mul_ui(s,s,4,prec);
    
    //左边I_c部分
    arb_sin(t,tau,prec);
    arb_neg(t,t);
    arb_mul(t,t,tau,prec);
    
    arb_mul(res,s,t,prec);
    
    arb_abs(w,res);
    if( arb_lt(w,INT_MIN_INTERVAL) ) //积分到正无穷时，会有些过小的结果
    {
        arb_zero(res);
        //arb_printn(res,25,0);printf("\n");
    }
    
    arb_get_mid_arb(res,res); //消误差
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(x_tau);
    arb_clear(y_tau);
    arb_clear(T_x_tau);
    arb_clear(T_y_tau);
    arb_clear(T_x_tau_prime);
    arb_clear(T_y_tau_prime);
    
    return 0;
}



static int interior_GW_I_s_func(arb_t res, const arb_t tau, void* param, const slong order, slong prec)
{
    arb_t s,t,w,x_tau,y_tau,T_x_tau,T_y_tau,T_x_tau_prime,T_y_tau_prime;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(x_tau);
    arb_init(y_tau);
    arb_init(T_x_tau);
    arb_init(T_y_tau);
    arb_init(T_x_tau_prime);
    arb_init(T_y_tau_prime);
    
    //接收传入参数
    struct Func_transfer_parameter *func_x_y; //这里不需要手动分配，只需将其指向传入的指针即可
    func_x_y=param;
    
    arb_mul(x_tau,func_x_y->p_1,tau, prec); //x*τ
    arb_mul(y_tau,func_x_y->p_2,tau, prec); //y*τ
    
    Radiation_transfer_function(T_x_tau,x_tau,NULL,0,prec); //T(x*τ)
    Radiation_transfer_function(T_y_tau,y_tau,NULL,0,prec); //T(y*τ)
    
    Radiation_transfer_function(T_x_tau_prime,x_tau,NULL,1,prec); //T'(x*τ)
    Radiation_transfer_function(T_y_tau_prime,y_tau,NULL,1,prec); //T'(y*τ)
    
    
    //积分函数形式为$ \tau( - \sin\tau)\cdot4\Big\{2\mathcal{T}(x\tau)\mathcal{T}(y\tau) + \Big[\mathcal{T}(x\tau) + x\tau \mathcal{T}'(x\tau)\Big]\Big[\mathcal{T}(y\tau) + y\tau \mathcal{T}'(y\tau)\Big]\Big\} $
    //右边 I_c 与 I_s 共同部分
    arb_mul(s,x_tau,T_x_tau_prime,prec);
    arb_add(s,s,T_x_tau,prec);
    
    arb_mul(t,y_tau,T_y_tau_prime,prec);
    arb_add(t,t,T_y_tau,prec);
    arb_mul(s,s,t,prec);
    
    arb_mul(t,T_x_tau,T_y_tau,prec);
    arb_mul_ui(t,t,2,prec);
    arb_add(s,s,t,prec);
    arb_mul_ui(s,s,4,prec);
    
    //左边I_s部分
    arb_cos(t,tau,prec);
    arb_mul(t,t,tau,prec);
    
    arb_mul(res,s,t,prec);
    
    arb_abs(w,res);
    if( arb_lt(w,INT_MIN_INTERVAL) ) //积分到正无穷时，会有些过小的结果
    {
        arb_zero(res);
        //arb_printn(res,25,0);printf("\n");
    }
    
    arb_get_mid_arb(res,res); //消误差
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(x_tau);
    arb_clear(y_tau);
    arb_clear(T_x_tau);
    arb_clear(T_y_tau);
    arb_clear(T_x_tau_prime);
    arb_clear(T_y_tau_prime);
    
    return 0;
}


//积分表达式
int GW_I_c_func(arb_t res, const arb_t x, const arb_t y, const arb_t k, slong prec)
{
    arb_t s,a,b,e;
    arb_init(s);
    arb_init(a);
    arb_init(b);
    arb_init(e);
    
    //积分上下限
    //arb_one(a); //积分下限，为 1/k or 0 or 1，这里我们设为1
    //arb_pos_inf(b); //积分上限 +∞
    //arb_set_str(e,"1E-25",prec);
    if( arb_is_zero(Int_GW_I_func_min) || arb_is_one(Int_GW_I_func_min) )
    {
        arb_set(a,Int_GW_I_func_min);
    }else
    {
        arb_inv(a,k,prec); //当积分下限不为0或1时，设为 1/k
    }
    
    
    //其中x和y是积分函数的参数，需传入
    //这里，对于结构体 Func_transfer_parameter 需手动分配内存
    struct Func_transfer_parameter *func_x_y = (struct Func_transfer_parameter *)calloc(1,sizeof(struct Func_transfer_parameter));
    
    arb_init(func_x_y->p_1);//使用arb_t变量前初始化
    arb_init(func_x_y->p_2);
    
    arb_set(func_x_y->p_1,x); //设定传递参数
    arb_set(func_x_y->p_2,y);
    
    //注意，这里因为积分上限是 +∞ 故而要使用 Double Exponential积分
    Integration_arb(s, interior_GW_I_c_func, func_x_y, 0,
                    a, Int_GW_I_func_max, Int_GW_I_func_precision,
                    Int_GW_I_func_iterate_min, Int_GW_I_func_iterate_max,
                    prec);
    
    arb_set(res,s);
    
    arb_clear(s);
    arb_clear(a);
    arb_clear(b);
    arb_clear(e);
    
    arb_clear(func_x_y->p_1);
    arb_clear(func_x_y->p_2);
    free(func_x_y); //手动释放自定义结构体内存
    
    return 0;
}


int GW_I_s_func(arb_t res, const arb_t x, const arb_t y, const arb_t k, slong prec)
{
    arb_t s,a,b,e;
    arb_init(s);
    arb_init(a);
    arb_init(b);
    arb_init(e);
    
    //积分上下限[a,b]
    //arb_one(a); //积分下限，为 1/k or 0 or 1，这里我们设为1
    //arb_pos_inf(b); //积分上限 +∞
    //arb_set_str(e,"1E-25",prec);
    if( arb_is_zero(Int_GW_I_func_min) || arb_is_one(Int_GW_I_func_min) )
    {
        arb_set(a,Int_GW_I_func_min);
    }else
    {
        arb_inv(a,k,prec); //当积分下限不为0或1时，设为 1/k
    }
    
    //其中x和y是积分函数的参数，需传入
    //这里，对于结构体 Func_transfer_parameter 需手动分配内存
    struct Func_transfer_parameter *func_x_y = (struct Func_transfer_parameter *)calloc(1,sizeof(struct Func_transfer_parameter));
    
    arb_init(func_x_y->p_1);//使用arb_t变量前初始化
    arb_init(func_x_y->p_2);
    
    arb_set(func_x_y->p_1,x); //设定传递参数
    arb_set(func_x_y->p_2,y);
    
    //注意，这里因为积分上限是 +∞ 故而要使用 Double Exponential积分
    Integration_arb(s, interior_GW_I_s_func, func_x_y, 0,
                    a, Int_GW_I_func_max, Int_GW_I_func_precision,
                    Int_GW_I_func_iterate_min, Int_GW_I_func_iterate_max,
                    prec);
    
    arb_set(res,s);
    
    arb_clear(s);
    arb_clear(a);
    arb_clear(b);
    arb_clear(e);
    
    arb_clear(func_x_y->p_1);
    arb_clear(func_x_y->p_2);
    free(func_x_y); //手动释放自定义结构体内存
    
    return 0;
}


//
//解析部分
//


//解析解函数的辅助函数
static void A_x_help_I(arb_t res, const arb_t x, const arb_t tau, slong prec)
{
    arb_t s,t,w;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //A(x)=1/τ^3 *[2+(x^2-1)*τ^2]
    arb_sqr(s,x,prec);
    arb_sub_ui(s,s,1,prec);
    arb_sqr(t,tau,prec);
    arb_mul(s,s,t,prec);
    arb_add_ui(s,s,2,prec);
    
    arb_pow_ui(t,tau,3,prec);
    arb_div(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
}

static void B_x_help_I(arb_t res, const arb_t x, const arb_t tau, slong prec)
{
    arb_t s,t,w;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //B(x)=1/τ^3 *[6+(x^2-1)*τ^2]
    arb_sqr(s,x,prec);
    arb_sub_ui(s,s,1,prec);
    arb_sqr(t,tau,prec);
    arb_mul(s,s,t,prec);
    arb_add_ui(s,s,6,prec);
    
    arb_pow_ui(t,tau,3,prec);
    arb_div(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
}

static void C_x_y_help_I(arb_t res, const arb_t x, const arb_t y, const arb_t tau, slong prec)
{
    arb_t s,t,w;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //C(x,y)=1/τ^4 *[6-(1+2*x^2-y^2)*τ^2]
    arb_sqr(s,x,prec);
    arb_mul_ui(s,s,2,prec);
    arb_add_ui(s,s,1,prec);
    arb_sqr(t,y,prec);
    arb_sub(s,s,t,prec);
    
    arb_sqr(t,tau,prec);
    arb_mul(s,s,t,prec);
    arb_neg(s,s);
    arb_add_ui(s,s,6,prec);
    
    arb_pow_ui(t,tau,4,prec);
    arb_div(res,s,t,prec);
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
}





//解析函数表达式
static void interior_GW_I_c_func_analyze(arb_t res, const arb_t d, const arb_t s, const arb_t tau, slong prec)
{
    arb_t z,t,w,q,d_2,s_2,tau_2,d_tau,s_tau,A_s,A_d,B_s,B_d,C_d_s,C_s_d;
    arb_t sin_tau,sin_d_tau,sin_s_tau,cos_tau,cos_d_tau,cos_s_tau;
    
    arb_init(z);
    arb_init(t);
    arb_init(w);
    arb_init(q);
    arb_init(d_2);
    arb_init(s_2);
    arb_init(tau_2);
    arb_init(d_tau);
    arb_init(s_tau);
    arb_init(A_s);
    arb_init(A_d);
    arb_init(B_s);
    arb_init(B_d);
    arb_init(C_d_s);
    arb_init(C_s_d);
    arb_init(sin_tau);
    arb_init(sin_d_tau);
    arb_init(sin_s_tau);
    arb_init(cos_tau);
    arb_init(cos_d_tau);
    arb_init(cos_s_tau);
    
    arb_sqr(d_2,d,prec); //d^2
    arb_sqr(s_2,s,prec); //s^2
    arb_sqr(tau_2,tau,prec); //τ^2
    
    arb_mul(d_tau,d,tau,prec); //d*τ
    arb_mul(s_tau,s,tau,prec); //s*τ
    
    
    if( !arb_is_finite(tau) ) //τ=+∞时
    {
        //-36*π* (s^2+d^2-2)^2 / (s^2-d^2)^3 * Θ(s-1)
        arb_const_pi(t,prec);
        arb_mul_ui(t,t,36,prec);
        arb_neg(t,t);
        
        arb_add(w,s_2,d_2,prec);
        arb_sub_ui(w,w,2,prec);
        arb_sqr(w,w,prec);
        arb_sub(z,s_2,d_2,prec);
        arb_pow_ui(z,z,3,prec);
        arb_div(w,w,z,prec);
        arb_mul(t,t,w,prec);
        
        arb_sub_ui(w,s,1,prec);
        Heaviside_Theta_function(z,w,prec);
        arb_mul(res,t,z,prec);
        
    }else if ( arb_is_zero(tau) ) //τ=0时
    {
        arb_zero(res);
    }else
    {
        //具体形式，参看1804.07732 中的 D.1
        A_x_help_I(A_s,s,tau,prec);//获得各函数值
        A_x_help_I(A_d,d,tau,prec);
        B_x_help_I(B_s,s,tau,prec);
        B_x_help_I(B_d,d,tau,prec);
        C_x_y_help_I(C_d_s,d,s,tau,prec);
        C_x_y_help_I(C_s_d,s,d,tau,prec);
        
        arb_sin(sin_tau,tau,prec); //使用值
        arb_sin(sin_d_tau,d_tau,prec);
        arb_sin(sin_s_tau,s_tau,prec);
        arb_cos(cos_tau,tau,prec);
        arb_cos(cos_d_tau,d_tau,prec);
        arb_cos(cos_s_tau,s_tau,prec);
        
        //大括号内，第一行左
        arb_mul(z,A_s,cos_tau,prec);
        arb_mul(w,C_d_s,sin_tau,prec);
        arb_add(z,z,w,prec);
        arb_mul(z,z,cos_d_tau,prec);
        
        //大括号内，第一行右
        arb_mul(t,A_d,cos_tau,prec);
        arb_mul(w,C_s_d,sin_tau,prec);
        arb_add(t,t,w,prec);
        arb_mul(t,t,cos_s_tau,prec);
        arb_sub(z,z,t,prec);
        
        //大括号内，第二行左
        arb_inv(t,tau_2,prec);
        arb_mul_ui(t,t,2,prec);
        arb_mul(t,t,cos_tau,prec);
        arb_mul(w,B_s,sin_tau,prec);
        arb_add(t,t,w,prec);
        arb_mul(t,t,d,prec);
        arb_mul(t,t,sin_d_tau,prec);
        arb_add(z,z,t,prec);
        
        //大括号内，第二行右
        arb_inv(t,tau_2,prec);
        arb_mul_ui(t,t,2,prec);
        arb_mul(t,t,cos_tau,prec);
        arb_mul(w,B_d,sin_tau,prec);
        arb_add(t,t,w,prec);
        arb_mul(t,t,s,prec);
        arb_mul(t,t,sin_s_tau,prec);
        arb_sub(z,z,t,prec);
        
        //大括号内，第三行左系数部分
        arb_add(t,s_2,d_2,prec);
        arb_sub_ui(t,t,2,prec);
        arb_sqr(t,t,prec);
        arb_div_ui(t,t,8,prec);
        
        //大括号内，第三行右边括号内
        arb_sub(q,tau,s_tau,prec);
        arb_hypgeom_si(q,q,prec); //Si(x)
        arb_add(w,tau,s_tau,prec);
        arb_hypgeom_si(w,w,prec);
        arb_add(q,q,w,prec);
        
        arb_add(w,tau,d_tau,prec);
        arb_hypgeom_si(w,w,prec);
        arb_sub(q,q,w,prec);
        
        arb_sub(w,tau,d_tau,prec);
        arb_hypgeom_si(w,w,prec);
        arb_sub(q,q,w,prec);
        
        arb_mul(t,t,q,prec);
        arb_add(z,z,t,prec);
        
        //大括号前面系数
        arb_sub(t,s_2,d_2,prec);
        arb_pow_ui(t,t,3,prec);
        arb_inv(t,t,prec);
        arb_mul_ui(t,t,288,prec);
        
        arb_mul(res,t,z,prec);
    }
    
    arb_clear(z);
    arb_clear(t);
    arb_clear(w);
    arb_clear(q);
    arb_clear(d_2);
    arb_clear(s_2);
    arb_clear(tau_2);
    arb_clear(d_tau);
    arb_clear(s_tau);
    arb_clear(A_s);
    arb_clear(A_d);
    arb_clear(B_s);
    arb_clear(B_d);
    arb_clear(C_d_s);
    arb_clear(C_s_d);
    arb_clear(sin_tau);
    arb_clear(sin_d_tau);
    arb_clear(sin_s_tau);
    arb_clear(cos_tau);
    arb_clear(cos_d_tau);
    arb_clear(cos_s_tau);
}


static void interior_GW_I_s_func_analyze(arb_t res, const arb_t d, const arb_t s, const arb_t tau, slong prec)
{
    arb_t z,t,w,q,d_2,s_2,tau_2,d_tau,s_tau,A_s,A_d,B_s,B_d,C_d_s,C_s_d;
    arb_t sin_tau,sin_d_tau,sin_s_tau,cos_tau,cos_d_tau,cos_s_tau;
    
    arb_init(z);
    arb_init(t);
    arb_init(w);
    arb_init(q);
    arb_init(d_2);
    arb_init(s_2);
    arb_init(tau_2);
    arb_init(d_tau);
    arb_init(s_tau);
    arb_init(A_s);
    arb_init(A_d);
    arb_init(B_s);
    arb_init(B_d);
    arb_init(C_d_s);
    arb_init(C_s_d);
    arb_init(sin_tau);
    arb_init(sin_d_tau);
    arb_init(sin_s_tau);
    arb_init(cos_tau);
    arb_init(cos_d_tau);
    arb_init(cos_s_tau);
    
    arb_sqr(d_2,d,prec); //d^2
    arb_sqr(s_2,s,prec); //s^2
    arb_sqr(tau_2,tau,prec); //τ^2
    
    arb_mul(d_tau,d,tau,prec); //d*τ
    arb_mul(s_tau,s,tau,prec); //s*τ
    
    
    if( !arb_is_finite(tau) ) //τ=+∞时
    {
        arb_zero(res);
        
    }else if ( arb_is_zero(tau) ) //τ=0时
    {
        //36* (s^2+d^2-2) / (s^2-d^2)^2 * [ (s^2+d^2-2)/(s^2-d^2)^2 * log( (1-d^2)/|s^2-1|) +2 ]
        
        arb_add(w,s_2,d_2,prec);//前面系数部分
        arb_sub_ui(w,w,2,prec);
        arb_sub(q,s_2,d_2,prec);
        arb_div(z,w,q,prec); //此结果后面要用
        arb_div(w,z,q,prec);
        arb_mul_ui(w,w,36,prec);
        
        arb_neg(q,d_2); //后面括号部分
        arb_add_ui(q,q,1,prec);
        arb_sub_ui(t,s_2,1,prec);
        arb_div(q,q,t,prec);
        arb_abs(q,q); //取对数前取绝对值
        arb_log(q,q,prec);
        arb_mul(q,q,z,prec); //使用前面的结果
        arb_add_ui(q,q,2,prec);
        
        arb_mul(res,w,q,prec);
    }else
    {
        //具体形式，参看1804.07732 中的 D.2
        A_x_help_I(A_s,s,tau,prec);//获得各函数值
        A_x_help_I(A_d,d,tau,prec);
        B_x_help_I(B_s,s,tau,prec);
        B_x_help_I(B_d,d,tau,prec);
        C_x_y_help_I(C_d_s,d,s,tau,prec);
        C_x_y_help_I(C_s_d,s,d,tau,prec);
        
        arb_sin(sin_tau,tau,prec); //使用值
        arb_sin(sin_d_tau,d_tau,prec);
        arb_sin(sin_s_tau,s_tau,prec);
        arb_cos(cos_tau,tau,prec);
        arb_cos(cos_d_tau,d_tau,prec);
        arb_cos(cos_s_tau,s_tau,prec);
        
        //大括号内，第一行左
        arb_mul(z,A_s,sin_tau,prec);
        arb_mul(w,C_d_s,cos_tau,prec);
        arb_sub(z,z,w,prec);
        arb_mul(z,z,cos_d_tau,prec);
        
        //大括号内，第一行右
        arb_mul(t,A_d,sin_tau,prec);
        arb_neg(t,t);
        arb_mul(w,C_s_d,cos_tau,prec);
        arb_add(t,t,w,prec);
        arb_mul(t,t,cos_s_tau,prec);
        arb_add(z,z,t,prec);
        
        //大括号内，第二行左
        arb_inv(t,tau_2,prec);
        arb_mul_ui(t,t,2,prec);
        arb_mul(t,t,sin_tau,prec);
        arb_mul(w,B_s,cos_tau,prec);
        arb_sub(t,t,w,prec);
        arb_mul(t,t,d,prec);
        arb_mul(t,t,sin_d_tau,prec);
        arb_add(z,z,t,prec);
        
        //大括号内，第二行右
        arb_inv(t,tau_2,prec);
        arb_mul_ui(t,t,2,prec);
        arb_neg(t,t);
        arb_mul(t,t,sin_tau,prec);
        arb_mul(w,B_d,cos_tau,prec);
        arb_add(t,t,w,prec);
        arb_mul(t,t,s,prec);
        arb_mul(t,t,sin_s_tau,prec);
        arb_add(z,z,t,prec);
        
        //大括号内，第三行左系数部分
        arb_add(t,s_2,d_2,prec);
        arb_sub_ui(t,t,2,prec);
        arb_sqr(t,t,prec);
        arb_div_ui(t,t,8,prec);
        
        //大括号内，第三行右边括号内
        arb_add(q,tau,d_tau,prec);
        arb_hypgeom_ci(q,q,prec); //Ci(x) 
        arb_sub(w,tau,d_tau,prec);
        arb_hypgeom_ci(w,w,prec);
        arb_add(q,q,w,prec);
        
        arb_sub(w,s_tau,tau,prec);
        arb_abs(w,w);
        arb_hypgeom_ci(w,w,prec);
        arb_sub(q,q,w,prec);
        
        arb_add(w,tau,s_tau,prec);
        arb_hypgeom_ci(w,w,prec);
        arb_sub(q,q,w,prec);
        
        arb_mul(t,t,q,prec);
        arb_add(z,z,t,prec);
        
        //大括号前面系数
        arb_sub(t,s_2,d_2,prec);
        arb_pow_ui(t,t,3,prec);
        arb_inv(t,t,prec);
        arb_mul_ui(t,t,288,prec);
        
        arb_mul(res,t,z,prec);
    }
    
    
    arb_clear(z);
    arb_clear(t);
    arb_clear(w);
    arb_clear(q);
    arb_clear(d_2);
    arb_clear(s_2);
    arb_clear(tau_2);
    arb_clear(d_tau);
    arb_clear(s_tau);
    arb_clear(A_s);
    arb_clear(A_d);
    arb_clear(B_s);
    arb_clear(B_d);
    arb_clear(C_d_s);
    arb_clear(C_s_d);
    arb_clear(sin_tau);
    arb_clear(sin_d_tau);
    arb_clear(sin_s_tau);
    arb_clear(cos_tau);
    arb_clear(cos_d_tau);
    arb_clear(cos_s_tau);
}


//解析解
int GW_I_c_func_analyze(arb_t res, const arb_t x, const arb_t y, const arb_t k, slong t_prec)
{
    arb_t q,w,a,b,d,s;
    arb_init(q);
    arb_init(w);
    arb_init(a);
    arb_init(b);
    arb_init(d);
    arb_init(s);
    
    //这里，临时的提升一下计算的精度
    slong prec;
    prec=t_prec/2*3; //提升 3/2 倍
    
    //积分上下限[a,b]
    //arb_one(a); //积分下限，为 1/k or 0 or 1，这里我们设为1
    //arb_pos_inf(b); //积分上限 +∞
    //arb_set_str(b,"3",prec);
    if( arb_is_zero(Int_GW_I_func_min) || arb_is_one(Int_GW_I_func_min) )
    {
        arb_set(a,Int_GW_I_func_min);
    }else
    {
        arb_inv(a,k,prec); //当积分下限不为0或1时，设为 1/k
    }
    
    
    //变量替换 (x,y)-->(d,s)
    arb_sub(q,x,y,prec); // 1/sqrt(3) * |x-y|
    arb_abs(q,q);
    arb_one(w);
    arb_mul_ui(w,w,3,prec);
    arb_sqrt(w,w,prec);
    arb_div(d,q,w,prec);
    
    arb_add(q,x,y,prec);
    arb_div(s,q,w,prec);
    
    //利用原函数求积分
    interior_GW_I_c_func_analyze(q,d,s,a,prec);
    interior_GW_I_c_func_analyze(w,d,s,Int_GW_I_func_max,prec);
    
    arb_sub(res,w,q,prec);
    
    arb_clear(q);
    arb_clear(w);
    arb_clear(a);
    arb_clear(b);
    arb_clear(d);
    arb_clear(s);
    
    return 0;
}


int GW_I_s_func_analyze(arb_t res, const arb_t x, const arb_t y, const arb_t k, slong t_prec)
{
    arb_t q,w,a,b,d,s;
    arb_init(q);
    arb_init(w);
    arb_init(a);
    arb_init(b);
    arb_init(d);
    arb_init(s);
    
    //这里，临时的提升一下计算的精度
    slong prec;
    prec=t_prec/2*3; //提升 3/2 倍
    
    //积分上下限[a,b]
    //arb_one(a); //积分下限，为 1/k or 0 or 1，这里我们设为1
    //arb_pos_inf(b); //积分上限 +∞
    //arb_set_str(b,"3",prec);
    if( arb_is_zero(Int_GW_I_func_min) || arb_is_one(Int_GW_I_func_min) )
    {
        arb_set(a,Int_GW_I_func_min);
    }else
    {
        arb_inv(a,k,prec); //当积分下限不为0或1时，设为 1/k
    }
    
    
    //变量替换 (x,y)-->(d,s)
    arb_sub(q,x,y,prec); // 1/sqrt(3) * |x-y|
    arb_abs(q,q);
    arb_one(w);
    arb_mul_ui(w,w,3,prec);
    arb_sqrt(w,w,prec);
    arb_div(d,q,w,prec);
    
    arb_add(q,x,y,prec);
    arb_div(s,q,w,prec);
    
    //利用原函数求积分
    interior_GW_I_s_func_analyze(q,d,s,a,prec);
    interior_GW_I_s_func_analyze(w,d,s,Int_GW_I_func_max,prec);
    
    arb_sub(res,w,q,prec);
    
    arb_clear(q);
    arb_clear(w);
    arb_clear(a);
    arb_clear(b);
    arb_clear(d);
    arb_clear(s);
    
    return 0;
}

