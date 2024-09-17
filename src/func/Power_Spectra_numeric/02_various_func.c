#include "Inflation_odes_slove.h"
#include <stdlib.h>  

//在此处，定义了微分方程计算所需函数


//求解背景扰动方程所需函数
void Inflation_Smoothing_Step_Function(arb_t res, const arb_t x, const arb_t smoothness, slong prec)
{
    arb_t s,t,w;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    arb_mul(t,x,smoothness,prec);
    arb_one(s);
    
    if( arb_is_negative(t) )
    {
        arb_zero(res);
    }else if( arb_gt(t,s) )
    {
        arb_one(res);
    }else
    {
        /*
         *       //(t**2) * (3 - 2.*t)
         *       arb_sqr(s,t,prec);
         *       arb_mul_ui(w,t,2,prec);
         *       arb_neg(w,w);
         *       arb_add_ui(w,w,3,prec);
         *       arb_mul(res,s,w,prec);
         */
        
        //(t**3) * (10. -15.*t + 6.*t**2.)
        arb_sqr(s,t,prec);
        arb_mul_ui(s,s,6,prec);
        arb_mul_si(w,t,-15,prec);
        arb_add(s,s,w,prec);
        arb_add_ui(s,s,10,prec);
        
        arb_pow_ui(w,t,3,prec);
        arb_mul(res,s,w,prec);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
}

void Inflation_S_prime(arb_t res, const arb_t x, const arb_t smoothness, slong prec)
{
    arb_t s,t,w;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    arb_mul(t,x,smoothness,prec);
    arb_one(s);
    
    if( arb_is_negative(t) )
    {
        arb_zero(res);
    }else if( arb_gt(t,s) )
    {
        arb_zero(res);
    }else
    {
        /*
         *        //6.*(t) * (1. - t)
         *        arb_mul_ui(s,t,6,prec);
         *        arb_neg(w,t);
         *        arb_add_ui(w,w,1,prec);
         *        arb_mul(res,s,w,prec);
         */
        
        //30.*(t**2.) * (1. - 2.*t +(t**2.))
        arb_sqr(w,t,prec);
        
        arb_mul_si(s,t,-2,prec);
        arb_add(s,s,w,prec);
        arb_add_si(s,s,1,prec);
        
        arb_mul_si(w,w,30,prec);
        arb_mul(res,w,s,prec);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
}

void Inflation_S_pp(arb_t res, const arb_t x, const arb_t smoothness, slong prec)
{
    arb_t s,t,w;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    arb_mul(t,x,smoothness,prec);
    arb_one(s);
    
    if( arb_is_negative(t) )
    {
        arb_zero(res);
    }else if( arb_gt(t,s) )
    {
        arb_zero(res);
    }else
    {
        //60.*(t) * (1. - 3.*t +2.*(t**2.))
        arb_sqr(w,t,prec);
        arb_mul_si(w,w,2,prec);
        arb_mul_si(s,t,-3,prec);
        arb_add(s,s,w,prec);
        arb_add_si(s,s,1,prec);
        
        arb_mul_si(w,t,60,prec);
        arb_mul(res,w,s,prec);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
}


void Inflation_V_phi(arb_t res, const arb_t phi, slong prec)
{
    arb_t s,t,w,phi_small,phi_large,smooth_small,smooth_large;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(phi_small);
    arb_init(phi_large);
    arb_init(smooth_small);
    arb_init(smooth_large);
    
    //phi after the step 
    arb_sub(phi_small,phi,Inf_Phi_e,prec);
    arb_sub(phi_large,phi,Inf_Phi_s,prec);
    
    arb_neg(w,phi_small);
    Inflation_Smoothing_Step_Function(smooth_small, w, Inf_Lambda, prec);
    Inflation_Smoothing_Step_Function(smooth_large, phi_large, Inf_Lambda, prec);
    
    //V = (V0 + (Delta_V + sqrt_E2 * phi_small + Eta_2 * phi_small**2) * smooth_small + 
    //     sqrt_E1 * phi_large * smooth_large)
    
    //中间部分
    //(Delta_V + sqrt_E2 * phi_small + Eta_2 * phi_small**2) * smooth_small
    arb_mul(s,Inf_sqrt_E2,phi_small,prec);
    arb_sqr(t,phi_small,prec);
    arb_mul(t,t,Inf_Eta_2,prec);
    arb_add(s,s,t,prec);
    arb_add(s,s,Inf_Delta_V,prec);
    arb_mul(s,s,smooth_small,prec);
    
    arb_add(s,s,Inf_V0,prec); //加上最前面的 V0
    
    //未尾部分
    //sqrt_E1 * phi_large * smooth_large
    
    arb_mul(t,Inf_sqrt_E1,phi_large,prec);
    arb_mul(t,t,smooth_large,prec);
    
    arb_add(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(phi_small);
    arb_clear(phi_large);
    arb_clear(smooth_small);
    arb_clear(smooth_large);
}


//The first derivertive of V(\phi) with respect to the phi
void Inflation_V_phi_p(arb_t res, const arb_t phi, slong prec)
{
    arb_t s,t,w,phi_small,phi_large,smooth_small,smooth_large;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(phi_small);
    arb_init(phi_large);
    arb_init(smooth_small);
    arb_init(smooth_large);
    
    //phi after the step
    arb_sub(phi_small,phi,Inf_Phi_e,prec);
    arb_sub(phi_large,phi,Inf_Phi_s,prec);
    
    arb_neg(w,phi_small);
    Inflation_Smoothing_Step_Function(smooth_small, w, Inf_Lambda, prec);
    Inflation_Smoothing_Step_Function(smooth_large, phi_large, Inf_Lambda, prec);
    
    //((Delta_V + sqrt_E2 * phi_small + Eta_2 * phi_small**2) * (-Lambda * S_prime(-phi_small, Lambda)) +
    //      (sqrt_E2 + Eta_2 * phi_small * 2.) * smooth_small + 
    //      sqrt_E1 * smooth_large + sqrt_E1 * phi_large * S_prime(phi_large, Lambda) * Lambda)
    
    //前面
    //(Delta_V + sqrt_E2 * phi_small + Eta_2 * phi_small**2) * (-Lambda * S_prime(-phi_small, Lambda))
    arb_mul(s,Inf_sqrt_E2,phi_small,prec);
    arb_sqr(t,phi_small,prec);
    arb_mul(t,t,Inf_Eta_2,prec);
    arb_add(s,s,t,prec);
    arb_add(s,s,Inf_Delta_V,prec);
    
    arb_neg(t,phi_small);
    Inflation_S_prime(w, t, Inf_Lambda, prec);
    arb_mul(w,w,Inf_Lambda,prec);
    arb_neg(w,w);
    arb_mul(s,s,w,prec);
    
    //中间
    //(sqrt_E2 + Eta_2 * phi_small * 2.) * smooth_small + sqrt_E1 * smooth_large
    arb_mul(t,Inf_Eta_2,phi_small,prec);
    arb_mul_ui(t,t,2,prec);
    arb_add(t,t,Inf_sqrt_E2,prec);
    arb_mul(t,t,smooth_small,prec);
    arb_add(s,s,t,prec);
    
    arb_mul(t,Inf_sqrt_E1,smooth_large,prec);
    arb_add(s,s,t,prec);
    
    //最后
    //sqrt_E1 * phi_large * S_prime(phi_large, Lambda) * Lambda
    
    arb_mul(t,Inf_sqrt_E1,phi_large,prec);
    Inflation_S_prime(w, phi_large, Inf_Lambda, prec);
    arb_mul(t,t,w,prec);
    arb_mul(t,t,Inf_Lambda,prec);
    
    arb_add(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(phi_small);
    arb_clear(phi_large);
    arb_clear(smooth_small);
    arb_clear(smooth_large);
}


void Inflation_V_phi_pp(arb_t res, const arb_t phi, slong prec)
{
    arb_t s,t,w,phi_small,phi_large,smooth_small,smooth_large;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(phi_small);
    arb_init(phi_large);
    arb_init(smooth_small);
    arb_init(smooth_large);
    
    //phi after the step
    arb_sub(phi_small,phi,Inf_Phi_e,prec);
    arb_sub(phi_large,phi,Inf_Phi_s,prec);
    
    arb_neg(w,phi_small);
    Inflation_Smoothing_Step_Function(smooth_small, w, Inf_Lambda, prec);
    Inflation_Smoothing_Step_Function(smooth_large, phi_large, Inf_Lambda, prec);
    
    // ((Delta_V + sqrt_E2 * phi_small + Eta_2 * phi_small**2) * ((Lambda**2) * S_pp(-phi_small, Lambda)) +
    //      (sqrt_E2 + Eta_2 * phi_small*2.) * (-Lambda * S_prime(-phi_small, Lambda)) +
    //      ( Eta_2 * 2.) * smooth_small + 
    //      (sqrt_E2 + Eta_2 * phi_small * 2.) * S_prime(-phi_small, Lambda)*(-Lambda ) +
    //      sqrt_E1 * S_prime(phi_large, Lambda)*(Lambda) + 
    //      sqrt_E1 * S_prime(phi_large, Lambda) * Lambda + 
    //      sqrt_E1 * phi_large * S_pp(phi_large, Lambda) * (Lambda**2.))
    
    //前半部分
    //(Delta_V + sqrt_E2 * phi_small + Eta_2 * phi_small**2) * ((Lambda**2) * S_pp(-phi_small, Lambda))
    arb_mul(s,Inf_sqrt_E2,phi_small,prec);
    arb_sqr(t,phi_small,prec);
    arb_mul(t,t,Inf_Eta_2,prec);
    arb_add(s,s,t,prec);
    arb_add(s,s,Inf_Delta_V,prec);
    
    arb_neg(t,phi_small);
    Inflation_S_pp(w, t, Inf_Lambda, prec);
    arb_sqr(t,Inf_Lambda,prec);
    arb_mul(t,t,w,prec);
    arb_mul(s,s,t,prec);
    
    //中间 (sqrt_E2 + Eta_2 * phi_small*2.) * (-Lambda * S_prime(-phi_small, Lambda))
    arb_neg(t,phi_small);
    Inflation_S_prime(w, t, Inf_Lambda, prec);
    arb_mul(w,w,Inf_Lambda,prec);
    arb_neg(w,w);
    
    arb_mul(t,Inf_Eta_2,phi_small,prec);
    arb_mul_si(t,t,2,prec);
    arb_add(t,t,Inf_sqrt_E2,prec);
    
    arb_mul(w,t,w,prec);
    arb_add(s,s,w,prec);
    
    //中间 ( Eta_2 * 2.) * smooth_small
    arb_mul_si(t,Inf_Eta_2,2,prec);
    arb_mul(t,t,smooth_small,prec);
    arb_add(s,s,t,prec);
    
    //中间 (sqrt_E2 + Eta_2 * phi_small * 2.) * S_prime(-phi_small, Lambda)*(-Lambda )
    //这里，与前面一样
    arb_add(s,s,w,prec);
    
    //中间 sqrt_E1 * S_prime(phi_large, Lambda)*(Lambda)
    Inflation_S_prime(w, phi_large, Inf_Lambda, prec);
    arb_mul(w,w,Inf_sqrt_E1,prec);
    arb_mul(w,w,Inf_Lambda,prec);
    
    arb_add(s,s,w,prec);
    
    //中间 sqrt_E1 * S_prime(phi_large, Lambda) * Lambda
    //这里，与前面一样
    arb_add(s,s,w,prec);
    
    //最后 sqrt_E1 * phi_large * S_pp(phi_large, Lambda) * (Lambda**2.)
    
    arb_mul(t,Inf_sqrt_E1,phi_large,prec);
    Inflation_S_pp(w, phi_large, Inf_Lambda, prec);
    arb_mul(t,t,w,prec);
    arb_sqr(w,Inf_Lambda,prec);
    arb_mul(t,t,w,prec);
    
    arb_add(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(phi_small);
    arb_clear(phi_large);
    arb_clear(smooth_small);
    arb_clear(smooth_large);
}


//有效质量
void Inflation_m_eff(arb_t res, const arb_t phi_dot, const arb_t phi_ddot,
                const arb_t H, const arb_t H_dot, const arb_t V_phi_pp, slong prec)
{
    arb_t s,t,w;    
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //V_phi_pp - (1/Mpl)**2 * (3.*phi_dot**2. + (2*phi_dot*phi_ddot)/H - (H_dot*(phi_dot**2.))/(H**2))
    
    //(3.*phi_dot**2. + (2*phi_dot*phi_ddot)/H - (H_dot*(phi_dot**2.))/(H**2))
    arb_sqr(s,phi_dot,prec);
    arb_mul_ui(s,s,3,prec);
    
    arb_mul(t,phi_dot,phi_ddot,prec);
    arb_mul_si(t,t,2,prec);
    arb_div(t,t,H,prec);
    arb_add(s,s,t,prec);
    
    arb_sqr(t,phi_dot,prec);
    arb_mul(t,t,H_dot,prec);
    arb_sqr(w,H,prec);
    arb_div(t,t,w,prec);
    arb_sub(s,s,t,prec);
    
    arb_inv(t,Inf_Mpl,prec);
    arb_sqr(t,t,prec);
    arb_mul(s,s,t,prec);
    
    arb_sub(res,V_phi_pp,s,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
}

