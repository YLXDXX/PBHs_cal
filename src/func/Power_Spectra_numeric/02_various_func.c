#include "Inflation_odes_slove.h"
#include <stdlib.h>  

//在此处，定义了微分方程计算所需函数


//求解背景扰动方程所需函数
void Inflation_Smoothing_Step_Function(arb_t res, const arb_t x, const arb_t smoothness, const int flag, slong prec)
{
    //这里添加 flag 主要是为了区别调用光滑函数的对像，从而处理没有USR的情况
    // flag = 1 为 SR_1 调用
    // flag = 2 为SR_2 调用
    arb_t s,t,w,yp_0;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(yp_0);
    
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
        //对于没有 USR 的情况，光滑函数使用 Hermite interpolation 的多项式
        if( arb_is_zero(Inf_Delta_Phi_usr) && flag == 1 )
        {
            //Hermite interpolation 函数及其一阶导数在端点处相等
            //s(t)=t+t*(t-1)*( 1-2*t+(t-1)*yp_0 )
            //s'(t)=1+(2*t-1)*(1-2*t+(t-1)*yp_0)+t*(t-1)*(yp_0-2)
            //s''(t)=2*(1-2*t+(t-1)*yp_0)+2*(2*t-1)*(yp_0-2)
            // SR1 为一段直线，斜率为 sqrt(2*ε_1)
            //其中 yp_0 的取值为 yp_0=sqrt(2*ε_1)/ΔV/(-λ), 能保证一阶导连续
            
            //yp_0=sqrt(2*ε_1)/ΔV/(-λ)
            arb_div(yp_0,Inf_sqrt_E1,Inf_Delta_V,prec);
            arb_div(yp_0,yp_0,smoothness,prec);
            arb_neg(yp_0,yp_0);
            
            //s(t)=t+t*(t-1)*( 1-2*t+(t-1)*yp_0 )
            
            arb_sub_ui(s,t,1,prec); //( 1-2*t+(t-1)*yp_0 )
            arb_mul(s,s,yp_0,prec);
            arb_add_ui(s,s,1,prec);
            arb_mul_ui(w,t,2,prec);
            arb_sub(s,s,w,prec);
            
            arb_mul(s,s,t,prec);
            arb_sub_ui(w,t,1,prec);
            arb_mul(s,s,w,prec);
            
            arb_add(res,s,t,prec);
            
        }else if ( arb_is_zero(Inf_Delta_Phi_usr) && flag == 2 ) //取消USR后， SR1 没有光滑函数参与
        {
            if( arb_is_positive(x) ) //传过来的是 与 ϕ-ϕ_s 的差值，这里相当于一个阶越函数
            {
                arb_one(res);
            }else
            {
                arb_zero(res);
            }
        }else
        {
            //当有 USR 阶段，采取 (t**3) * (10. -15.*t + 6.*t**2.) 形式的光滑函数
            
            arb_sqr(s,t,prec);
            arb_mul_ui(s,s,6,prec);
            arb_mul_si(w,t,-15,prec);
            arb_add(s,s,w,prec);
            arb_add_ui(s,s,10,prec);
            
            arb_pow_ui(w,t,3,prec);
            arb_mul(res,s,w,prec);
        }
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(yp_0);
}

void Inflation_S_prime(arb_t res, const arb_t x, const arb_t smoothness, const int flag, slong prec)
{
    arb_t s,t,w,yp_0;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(yp_0);
    
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
        //对于没有 USR 的情况，光滑函数使用 Hermite interpolation 的多项式
        if( arb_is_zero(Inf_Delta_Phi_usr) && flag == 1 )
        {
            //yp_0=sqrt(2*ε_1)/ΔV/(-λ)
            arb_div(yp_0,Inf_sqrt_E1,Inf_Delta_V,prec);
            arb_div(yp_0,yp_0,smoothness,prec);
            arb_neg(yp_0,yp_0);
            
            //s'(t)=1+(2*t-1)*(1-2*t+(t-1)*yp_0)+t*(t-1)*(yp_0-2)
            
            arb_sub_ui(s,t,1,prec); //(1-2*t+(t-1)*yp_0)
            arb_mul(s,s,yp_0,prec);
            arb_add_ui(s,s,1,prec);
            arb_mul_ui(w,t,2,prec);
            arb_sub(s,s,w,prec);
            
            arb_sub_ui(w,w,1,prec);
            arb_mul(s,s,w,prec);
            arb_add_ui(s,s,1,prec);
            
            arb_sub_ui(w,t,1,prec); //t*(t-1)*(yp_0-2)
            arb_mul(w,w,t,prec);
            arb_sub_ui(yp_0,yp_0,2,prec); // yp_0 用掉
            arb_mul(w,w,yp_0,prec);
            
            arb_add(res,s,w,prec);
            
        }else if( arb_is_zero(Inf_Delta_Phi_usr) && flag == 2 ) //取消USR后， SR1 没有光滑函数参与
        {
            arb_zero(res);
        }else
        {
            //当有 USR 阶段，采取 (t**3) * (10. -15.*t + 6.*t**2.) 形式的光滑函数
            //一阶导 30.*(t**2.) * (1. - 2.*t +(t**2.))
            arb_sqr(w,t,prec);
            
            arb_mul_si(s,t,-2,prec);
            arb_add(s,s,w,prec);
            arb_add_si(s,s,1,prec);
            
            arb_mul_si(w,w,30,prec);
            arb_mul(res,w,s,prec);
        }
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(yp_0);
}

void Inflation_S_pp(arb_t res, const arb_t x, const arb_t smoothness, const int flag, slong prec)
{
    arb_t s,t,w,yp_0;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(yp_0);
    
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
        //对于没有 USR 的情况，光滑函数使用 Hermite interpolation 的多项式
        if( arb_is_zero(Inf_Delta_Phi_usr) && flag == 1 )
        {
            //yp_0=sqrt(2*ε_1)/ΔV/(-λ)
            arb_div(yp_0,Inf_sqrt_E1,Inf_Delta_V,prec);
            arb_div(yp_0,yp_0,smoothness,prec);
            arb_neg(yp_0,yp_0);
            
            //s''(t)=2*(1-2*t+(t-1)*yp_0)+2*(2*t-1)*(yp_0-2)
            
            arb_sub_ui(s,t,1,prec); //(1-2*t+(t-1)*yp_0)
            arb_mul(s,s,yp_0,prec);
            arb_add_ui(s,s,1,prec);
            arb_mul_ui(w,t,2,prec);
            arb_sub(s,s,w,prec);
            
            arb_mul_ui(s,s,2,prec);
            
            arb_sub_ui(w,w,1,prec); //2*(2*t-1)*(yp_0-2)
            arb_sub_ui(yp_0,yp_0,2,prec); // yp_0 用掉
            arb_mul(w,w,yp_0,prec);
            arb_mul_ui(w,w,2,prec);
            
            arb_add(res,s,w,prec);
            
        }else if( arb_is_zero(Inf_Delta_Phi_usr) && flag == 2 ) //取消USR后， SR1 没有光滑函数参与
        {
            arb_zero(res);
        }else
        {
            //当有 USR 阶段，采取 (t**3) * (10. -15.*t + 6.*t**2.) 形式的光滑函数
            //二阶导 60.*(t) * (1. - 3.*t +2.*(t**2.))
            arb_sqr(w,t,prec);
            arb_mul_si(w,w,2,prec);
            arb_mul_si(s,t,-3,prec);
            arb_add(s,s,w,prec);
            arb_add_si(s,s,1,prec);
            
            arb_mul_si(w,t,60,prec);
            arb_mul(res,w,s,prec);
        }
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(yp_0);
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
    Inflation_Smoothing_Step_Function(smooth_small, w, Inf_Lambda, 1, prec);
    Inflation_Smoothing_Step_Function(smooth_large, phi_large, Inf_Lambda, 2, prec);
    
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
    Inflation_Smoothing_Step_Function(smooth_small, w, Inf_Lambda, 1, prec);
    Inflation_Smoothing_Step_Function(smooth_large, phi_large, Inf_Lambda, 2, prec);
    
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
    Inflation_S_prime(w, t, Inf_Lambda, 1, prec);
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
    Inflation_S_prime(w, phi_large, Inf_Lambda, 2, prec);
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
    Inflation_Smoothing_Step_Function(smooth_small, w, Inf_Lambda, 1, prec);
    Inflation_Smoothing_Step_Function(smooth_large, phi_large, Inf_Lambda, 2, prec);
    
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
    Inflation_S_pp(w, t, Inf_Lambda, 1, prec);
    arb_sqr(t,Inf_Lambda,prec);
    arb_mul(t,t,w,prec);
    arb_mul(s,s,t,prec);
    
    //中间 (sqrt_E2 + Eta_2 * phi_small*2.) * (-Lambda * S_prime(-phi_small, Lambda))
    arb_neg(t,phi_small);
    Inflation_S_prime(w, t, Inf_Lambda, 1, prec);
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
    Inflation_S_prime(w, phi_large, Inf_Lambda, 2, prec);
    arb_mul(w,w,Inf_sqrt_E1,prec);
    arb_mul(w,w,Inf_Lambda,prec);
    
    arb_add(s,s,w,prec);
    
    //中间 sqrt_E1 * S_prime(phi_large, Lambda) * Lambda
    //这里，与前面一样
    arb_add(s,s,w,prec);
    
    //最后 sqrt_E1 * phi_large * S_pp(phi_large, Lambda) * (Lambda**2.)
    
    arb_mul(t,Inf_sqrt_E1,phi_large,prec);
    Inflation_S_pp(w, phi_large, Inf_Lambda, 2, prec);
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

