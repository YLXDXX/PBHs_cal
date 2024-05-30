#include "02_power_spectra.h"
#include <stdlib.h>
#include <arb_hypgeom.h>

//未取极限，未取平均
int I_func_v_u_x(arb_t res, const arb_t v, const arb_t u, const arb_t x, slong prec)
{
    arb_t s,t,w,q;
    arb_t v_u,u_x,v_x,v_add_u,v_sub_u,sqrt_3;
    arb_t ux_div_sqrt_3,vx_div_sqrt_3,v_2_add_u_2_sub_3;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(q);
    arb_init(v_u);
    arb_init(u_x);
    arb_init(v_x);
    arb_init(v_add_u);
    arb_init(v_sub_u);
    arb_init(sqrt_3);
    arb_init(ux_div_sqrt_3);
    arb_init(vx_div_sqrt_3);
    arb_init(v_2_add_u_2_sub_3);
    
    arb_mul(v_u,v,u,prec); //v*u
    arb_mul(u_x,u,x,prec); //u*x
    arb_mul(v_x,v,x,prec); //v*x
    arb_add(v_add_u,v,u,prec); //v+u
    arb_sub(v_sub_u,v,u,prec); //v-u
    arb_one(sqrt_3); //sqrt(3)
    arb_mul_ui(sqrt_3,sqrt_3,3,prec);
    arb_sqrt(sqrt_3,sqrt_3,prec);
    arb_div(ux_div_sqrt_3,u_x,sqrt_3,prec); // u*x/sqrt(3)
    arb_div(vx_div_sqrt_3,v_x,sqrt_3,prec); // v*x/sqrt(3)
    arb_sqr(t,u,prec); // v^2+u^2-3
    arb_sqr(w,v,prec);
    arb_add(v_2_add_u_2_sub_3,t,w,prec);
    arb_sub_ui(v_2_add_u_2_sub_3,v_2_add_u_2_sub_3,3,prec); 
    
    //函数具体形式，参见 1804.08577 (22)
    arb_set_str(s,"1E5",prec);
    if ( arb_gt(x,s) ) // x --> +∞时其渐近行为
    {
        //右边大括号中的左边
        arb_mul_ui(s,v_u,4,prec);
        arb_neg(s,s);
        
        arb_sqr(t,v_add_u,prec);
        arb_neg(t,t);
        arb_add_ui(t,t,3,prec);
        arb_sqr(w,v_sub_u,prec);
        arb_neg(w,w);
        arb_add_ui(w,w,3,prec);
        arb_div(t,t,w,prec);
        arb_abs(t,t);
        arb_log(t,t,prec);
        arb_mul(t,t,v_2_add_u_2_sub_3,prec);
        arb_add(s,s,t,prec);
        
        arb_sin(t,x,prec);
        arb_mul(s,s,t,prec);
        
        //右边大括号中的右边
        arb_sub(w,v_add_u,sqrt_3,prec);
        Heaviside_Theta_function(t,w,prec);
        
        arb_const_pi(w,prec);
        arb_mul(w,w,v_2_add_u_2_sub_3,prec);
        arb_mul(t,t,w,prec);
        arb_cos(w,w,prec);
        arb_mul(t,t,w,prec);
        
        arb_sub(s,s,t,prec); //两项相减
        
        //前面系数部分
        arb_mul_ui(t,v_2_add_u_2_sub_3,3,prec);
        arb_pow_ui(w,v_u,3,prec);
        arb_mul(w,w,x,prec);
        arb_mul_ui(w,w,4,prec);
        arb_div(t,t,w,prec);
        
        arb_mul(res,s,t,prec);
        
    }else if ( arb_is_zero(x) ) //x=0时
    {
        arb_zero(res);
    }else
    {
        //特大括号内容
        //第一个括号部分
        arb_mul(w,u_x,v_x,prec); //第一行括号左边
        arb_mul(s,w,v_2_add_u_2_sub_3,prec);
        arb_mul(s,s,x,prec);
        arb_sin(t,x,prec);
        arb_mul(s,s,t,prec);
        
        arb_mul_ui(t,w,6,prec); //第一行括号右边
        arb_cos(w,ux_div_sqrt_3,prec);
        arb_mul(t,t,w,prec);
        arb_cos(w,vx_div_sqrt_3,prec);
        arb_mul(t,t,w,prec);
        arb_sub(s,s,t,prec);
        
        arb_mul_ui(t,sqrt_3,6,prec); //第二行括号左边
        arb_mul(t,t,u_x,prec);
        arb_cos(w,ux_div_sqrt_3,prec);
        arb_sin(q,vx_div_sqrt_3,prec);
        arb_mul(w,w,q,prec);
        arb_mul(t,t,w,prec);
        arb_add(s,s,t,prec);
        
        arb_mul_ui(t,sqrt_3,6,prec); //第二行括号中间
        arb_mul(t,t,v_x,prec);
        arb_sin(w,ux_div_sqrt_3,prec);
        arb_cos(q,vx_div_sqrt_3,prec);
        arb_mul(w,w,q,prec);
        arb_mul(t,t,w,prec);
        arb_add(s,s,t,prec);
        
        arb_sqr(w,x,prec); //第二行括号右边
        arb_mul(t,v_2_add_u_2_sub_3,w,prec);
        arb_add_ui(t,t,6,prec);
        arb_mul_ui(t,t,3,prec);
        arb_sin(w,ux_div_sqrt_3,prec);
        arb_sin(q,vx_div_sqrt_3,prec);
        arb_mul(w,w,q,prec);
        arb_mul(t,t,w,prec);
        arb_sub(s,s,t,prec);
        
        arb_pow_ui(t,x,3,prec); //第一个括号前面系数部分
        arb_inv(t,t,prec);
        arb_mul_ui(t,t,4,prec);
        arb_neg(t,t);
        arb_mul(s,s,t,prec);
        
        //第二个括号部分中的 sin(x)*(...)
        arb_sub(t,vx_div_sqrt_3,ux_div_sqrt_3,prec); //第三行
        arb_neg(t,t);
        arb_add(t,t,x,prec);
        arb_hypgeom_ci(t,t,prec); //Ci(x)
        
        arb_sub(w,vx_div_sqrt_3,ux_div_sqrt_3,prec);
        arb_add(w,w,x,prec);
        arb_hypgeom_ci(w,w,prec); //Ci(x)
        arb_add(t,t,w,prec);
        
        arb_add(w,vx_div_sqrt_3,ux_div_sqrt_3,prec); //第四行
        arb_neg(w,w);
        arb_add(w,w,x,prec);
        arb_abs(w,w);
        arb_hypgeom_ci(w,w,prec); //Ci(x)
        arb_sub(t,t,w,prec);
        
        arb_add(w,vx_div_sqrt_3,ux_div_sqrt_3,prec);
        arb_add(w,w,x,prec);
        arb_hypgeom_ci(w,w,prec); //Ci(x)
        arb_sub(t,t,w,prec);
        
        arb_sqr(w,v_add_u,prec);
        arb_neg(w,w);
        arb_add_ui(w,w,3,prec);
        arb_sqr(q,v_sub_u,prec);
        arb_neg(q,q);
        arb_add_ui(q,q,3,prec);
        arb_div(w,w,q,prec);
        arb_abs(w,w);
        arb_log(w,w,prec);
        arb_add(t,t,w,prec);
        
        arb_sin(w,x,prec); //括号前面系数
        arb_mul(t,t,w,prec);
        
        //第二个括号部分中的 cos(x)*(...)
        arb_sub(w,vx_div_sqrt_3,ux_div_sqrt_3,prec); //第五行
        arb_neg(w,w);
        arb_add(w,w,x,prec);
        arb_hypgeom_si(w,w,prec); //Si(x)
        arb_neg(w,w);
        
        arb_sub(q,vx_div_sqrt_3,ux_div_sqrt_3,prec);
        arb_add(q,q,x,prec);
        arb_hypgeom_si(q,q,prec); //Si(x)
        arb_sub(w,w,q,prec);
        
        arb_add(q,vx_div_sqrt_3,ux_div_sqrt_3,prec); //第六行
        arb_neg(q,q);
        arb_add(q,q,x,prec);
        arb_hypgeom_si(q,q,prec); //Si(x)
        arb_add(w,w,q,prec);
        
        arb_add(q,vx_div_sqrt_3,ux_div_sqrt_3,prec);
        arb_add(q,q,x,prec);
        arb_hypgeom_si(q,q,prec); //Si(x)
        arb_add(w,w,q,prec);
        
        arb_cos(q,x,prec); //括号前面系数
        arb_mul(w,w,q,prec);
        
        arb_add(t,t,w,prec); //括号中两项相加
        
        arb_sqr(w,v_2_add_u_2_sub_3,prec); //括号前面系数
        arb_mul(t,t,w,prec);
        
        arb_add(s,s,t,prec); //括号中两项相加
        
        arb_pow_ui(w,v_u,3,prec); //大括号前面的系数
        arb_mul(w,w,x,prec);
        arb_mul_ui(w,w,4,prec);
        arb_inv(w,w,prec);
        arb_mul_ui(w,w,3,prec);
        
        arb_mul(res,s,w,prec);
    }
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(q);
    arb_clear(v_u);
    arb_clear(u_x);
    arb_clear(v_x);
    arb_clear(v_add_u);
    arb_clear(v_sub_u);
    arb_clear(sqrt_3);
    arb_clear(ux_div_sqrt_3);
    arb_clear(vx_div_sqrt_3);
    arb_clear(v_2_add_u_2_sub_3);
    
    return 0;
}




//取 x --> +∞ 且平方后取平均
int I_func_square_v_u_x_average(arb_t res, const arb_t v, const arb_t u, const arb_t x, slong prec)
{
    arb_t s,t,w,q;
    arb_t v_u,u_x,v_x,v_add_u,v_sub_u,sqrt_3;
    arb_t ux_div_sqrt_3,vx_div_sqrt_3,v_2_add_u_2_sub_3;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(q);
    arb_init(v_u);
    arb_init(u_x);
    arb_init(v_x);
    arb_init(v_add_u);
    arb_init(v_sub_u);
    arb_init(sqrt_3);
    arb_init(ux_div_sqrt_3);
    arb_init(vx_div_sqrt_3);
    arb_init(v_2_add_u_2_sub_3);
    
    arb_mul(v_u,v,u,prec); //v*u
    arb_mul(u_x,u,x,prec); //u*x
    arb_mul(v_x,v,x,prec); //v*x
    arb_add(v_add_u,v,u,prec); //v+u
    arb_sub(v_sub_u,v,u,prec); //v-u
    arb_one(sqrt_3); //sqrt(3)
    arb_mul_ui(sqrt_3,sqrt_3,3,prec);
    arb_sqrt(sqrt_3,sqrt_3,prec);
    arb_div(ux_div_sqrt_3,u_x,sqrt_3,prec); // u*x/sqrt(3)
    arb_div(vx_div_sqrt_3,v_x,sqrt_3,prec); // v*x/sqrt(3)
    arb_sqr(t,u,prec); // v^2+u^2-3
    arb_sqr(w,v,prec);
    arb_add(v_2_add_u_2_sub_3,t,w,prec);
    arb_sub_ui(v_2_add_u_2_sub_3,v_2_add_u_2_sub_3,3,prec); 
    
    //函数具体形式，参见 1804.08577 (26)
    arb_set_str(s,"1E5",prec);
    if ( arb_gt(x,s) ) // x --> +∞时其渐近行为
    {
        //右边大括号中的左边
        arb_mul_ui(s,v_u,4,prec);
        arb_neg(s,s);
        
        arb_sqr(t,v_add_u,prec);
        arb_neg(t,t);
        arb_add_ui(t,t,3,prec);
        arb_sqr(w,v_sub_u,prec);
        arb_neg(w,w);
        arb_add_ui(w,w,3,prec);
        arb_div(t,t,w,prec);
        arb_abs(t,t);
        arb_log(t,t,prec);
        arb_mul(t,t,v_2_add_u_2_sub_3,prec);
        arb_add(s,s,t,prec);
        
        arb_sqr(s,s,prec);
        
        //右边大括号中的右边
        arb_sub(w,v_add_u,sqrt_3,prec);
        Heaviside_Theta_function(t,w,prec);
        
        arb_const_pi(w,prec);
        arb_mul(w,w,v_2_add_u_2_sub_3,prec);
        arb_mul(t,t,w,prec);
        arb_sqr(t,t,prec);
        
        arb_add(s,s,t,prec); //两项相加
        
        //前面系数部分
        arb_mul_ui(t,v_2_add_u_2_sub_3,3,prec);
        arb_pow_ui(w,v_u,3,prec);
        arb_mul(w,w,x,prec);
        arb_mul_ui(w,w,4,prec);
        arb_div(t,t,w,prec);
        arb_sqr(t,t,prec);
        arb_div_ui(t,t,2,prec);
        
        arb_mul(res,s,t,prec);
        
    }else if ( arb_is_zero(x) ) //x=0时
    {
        arb_zero(res);
    }else 
    {
        printf("暂未实现\n");
        exit(1);
    }
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(q);
    arb_clear(v_u);
    arb_clear(u_x);
    arb_clear(v_x);
    arb_clear(v_add_u);
    arb_clear(v_sub_u);
    arb_clear(sqrt_3);
    arb_clear(ux_div_sqrt_3);
    arb_clear(vx_div_sqrt_3);
    arb_clear(v_2_add_u_2_sub_3);
    
    return 0;
}

