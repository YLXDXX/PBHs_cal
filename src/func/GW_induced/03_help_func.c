#include "03_help_func.h"
#include <stdlib.h>
#include <arb_hypgeom.h>

static void I_A_u_v(arb_t res, const arb_t u, const arb_t v, const arb_t u_2_v_2_sub_3, const arb_t uv, slong prec)
{
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    arb_mul_ui(s,u_2_v_2_sub_3,3,prec);
    
    arb_pow_ui(t,uv,3,prec);
    arb_mul_ui(t,t,4,prec);
    
    arb_div(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
}

static void I_B_u_v(arb_t res, const arb_t u, const arb_t v,const arb_t u_2_v_2_sub_3, const arb_t uv, slong prec)
{
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    arb_add(s,u,v,prec);
    arb_sqr(s,s,prec);
    arb_neg(s,s);
    arb_add_ui(s,s,3,prec);
    
    arb_sub(t,u,v,prec);
    arb_sqr(t,t,prec);
    arb_neg(t,t);
    arb_add_ui(t,t,3,prec);
    
    arb_div(s,s,t,prec);
    arb_abs(s,s);
    arb_log(s,s,prec);
    arb_mul(s,s,u_2_v_2_sub_3,prec);
    
    arb_mul_si(t,uv,-4,prec);
    
    arb_add(res,t,s,prec);
    
    arb_clear(s);
    arb_clear(t);
}



static void I_C_u_v(arb_t res, const arb_t u, const arb_t v, const arb_t u_2_v_2_sub_3, const arb_t uv, const arb_t sqrt_3, slong prec)
{
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    arb_add(s,u,v,prec);
    arb_sub(s,s,sqrt_3,prec);
    Heaviside_Theta_function(t,s,prec);
    
    arb_mul(res,u_2_v_2_sub_3,t,prec);
    
    arb_clear(s);
    arb_clear(t);
}


//u_1=u_2 and v_1=v_2 的简情况
static void I_RD_2_u_v_limited(arb_t res, const arb_t u, const arb_t v, slong prec) //取极限后，不含x
{
    arb_t s,t,I_A,I_B,I_C;
    arb_init(s);
    arb_init(t);
    arb_init(I_A);
    arb_init(I_B);
    arb_init(I_C);
    
    //此函数的具休形式见 2305.19950 (3.17)
    arb_t u_2_v_2_sub_3,uv,sqrt_3;
    arb_init(u_2_v_2_sub_3);
    arb_init(uv);
    arb_init(sqrt_3);
    
    //函数用值
    arb_one(sqrt_3);
    arb_mul_ui(sqrt_3,sqrt_3,3,prec);
    arb_sqrt(sqrt_3,sqrt_3,prec);
    
    arb_sqr(s,u,prec);//u^2+v^2-3
    arb_sqr(t,v,prec);
    arb_add(s,s,t,prec);
    arb_sub_ui(u_2_v_2_sub_3,s,3,prec);
    
    arb_mul(uv,u,v,prec); //u*v
    
    I_A_u_v(I_A,u,v,u_2_v_2_sub_3,uv,prec);
    arb_sqr(I_A,I_A,prec);
    
    I_B_u_v(I_B,u,v,u_2_v_2_sub_3,uv,prec);
    arb_sqr(I_B,I_B,prec);
    
    I_C_u_v(I_C,u,v,u_2_v_2_sub_3,uv,sqrt_3,prec);
    arb_sqr(I_C,I_C,prec);
    
    arb_sqr(s,Pi,prec);
    arb_mul(s,s,I_C,prec);
    arb_add(s,s,I_B,prec);
    
    arb_mul(s,s,I_A,prec);
    arb_div_ui(res,s,2,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(I_A);
    arb_clear(I_B);
    arb_clear(I_C);
    arb_clear(u_2_v_2_sub_3);
    arb_clear(uv);
    arb_clear(sqrt_3);
}

void I_RD_u_1_v_1_u_2_v_2_limited(arb_t res, const arb_t u_1, const arb_t v_1, const arb_t u_2, const arb_t v_2, slong prec) //取极限后，不含x
{
    arb_t s,t,I_A_1,I_B_1,I_C_1,I_A_2,I_B_2,I_C_2;
    arb_init(s);
    arb_init(t);
    arb_init(I_A_1);
    arb_init(I_B_1);
    arb_init(I_C_1);
    arb_init(I_A_2);
    arb_init(I_B_2);
    arb_init(I_C_2);
    
    //此函数的具休形式见 2305.19950 (3.17)
    arb_t u_2_v_2_sub_3,uv,sqrt_3;
    arb_init(u_2_v_2_sub_3);
    arb_init(uv);
    arb_init(sqrt_3);
    
    //函数用值
    arb_one(sqrt_3);
    arb_mul_ui(sqrt_3,sqrt_3,3,prec);
    arb_sqrt(sqrt_3,sqrt_3,prec);
    
    //u_1 v_1相关
    arb_sqr(s,u_1,prec);//u^2+v^2-3
    arb_sqr(t,v_1,prec);
    arb_add(s,s,t,prec);
    arb_sub_ui(u_2_v_2_sub_3,s,3,prec);
    
    arb_mul(uv,u_1,v_1,prec); //u*v
    
    I_A_u_v(I_A_1,u_1,v_1,u_2_v_2_sub_3,uv,prec);
    I_B_u_v(I_B_1,u_1,v_1,u_2_v_2_sub_3,uv,prec);
    I_C_u_v(I_C_1,u_1,v_1,u_2_v_2_sub_3,uv,sqrt_3,prec);

    //u_2 v_2相关
    arb_sqr(s,u_2,prec);//u^2+v^2-3
    arb_sqr(t,v_2,prec);
    arb_add(s,s,t,prec);
    arb_sub_ui(u_2_v_2_sub_3,s,3,prec);
    
    arb_mul(uv,u_2,v_2,prec); //u*v
    
    I_A_u_v(I_A_2,u_2,v_2,u_2_v_2_sub_3,uv,prec);
    I_B_u_v(I_B_2,u_2,v_2,u_2_v_2_sub_3,uv,prec);
    I_C_u_v(I_C_2,u_2,v_2,u_2_v_2_sub_3,uv,sqrt_3,prec);
    
    
    arb_sqr(s,Pi,prec);
    arb_mul(s,s,I_C_1,prec);
    arb_mul(s,s,I_C_2,prec);
    arb_mul(t,I_B_1,I_B_2,prec);
    arb_add(s,s,t,prec);
    
    arb_mul(s,s,I_A_1,prec);
    arb_mul(s,s,I_A_2,prec);
    arb_div_ui(res,s,2,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(I_A_1);
    arb_clear(I_B_1);
    arb_clear(I_C_1);
    arb_clear(I_A_2);
    arb_clear(I_B_2);
    arb_clear(I_C_2);
    arb_clear(u_2_v_2_sub_3);
    arb_clear(uv);
    arb_clear(sqrt_3);
}


void J_2_u_1_v_1_limited(arb_t res, const arb_t u_1, const arb_t v_1, slong prec) //取极限后，不含x
{
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    //此函数的具休形式见 2305.19950 (4.19)
    arb_add(s,v_1,u_1,prec);
    arb_sqr(s,s,prec);
    arb_sub_ui(s,s,1,prec);
    arb_sqr(s,s,prec);
    
    arb_sub(t,v_1,u_1,prec);
    arb_sqr(t,t,prec);
    arb_neg(t,t);
    arb_add_ui(t,t,1,prec);
    arb_sqr(t,t,prec);
    
    arb_mul(s,s,t,prec);
    arb_div_ui(s,s,64,prec);
    
    I_RD_2_u_v_limited(t,u_1,v_1,prec);
    arb_mul(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
}

void J_u_1_v_1_and_u_2_v_2_limited(arb_t res, const arb_t u_1, const arb_t v_1, const arb_t u_2, const arb_t v_2, slong prec) //取极限后，不含x
{
    
}

