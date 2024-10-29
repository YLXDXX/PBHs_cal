#include "Inflation_odes_slove.h"
#include <stdlib.h>  

//在此处，定义了与暴胀子势能相关的函数 
//包括：势能 V(ϕ)、分段势能的光滑连接函数 S(x;λ)


//分段势能光滑连接函数，及其导数
// S(x,λ) 其中 λ 用于控制光滑程度
static void Inflation_Potential_Smoothing_Function(arb_t res, const arb_t x, const arb_t smoothness,
                                            const slong order, const slong flag, slong prec)
{
    //这里添加 flag 主要是为了区别调用光滑函数的对像，从而处理没有USR的情况
    // flag = 1 为 phi_small 调用
    // flag = 2 为 phi_large 调用
    
    arb_t s,t,w,yp_0;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(yp_0);
    
    switch ( order ) //这里的 order 用于识别导数
    {
        case 0: //原函数
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
            
            break;
        case 1: //一阶导
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
            
            break;
        case 2: //二阶导
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
            
            break;
        default:
            printf("src/func/Power_Spectra_numeric/02_various_func.c Inflation_Potential_Smoothing_Function 阶数输入有误\n");
            exit(-1);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(yp_0);
}



//暴胀子的势能，及其导数 V(ϕ)
//对于多段的势能，分开定义各段的势能及其导数
//这里，考虑到光滑函数，使用有量纲的慢滚参数较为简单
//并且，势能的基点 V_0 并不写入每个分段函数
//在总势能那里统一处理，可简化与光滑函数相关的表达式

//第一段势能，及其导数
//V_0 在总势能那里统一处理
//第一个 SR 阶段，V_I(ϕ)=[ sqrt_E1*(ϕ-ϕ_s) + η_1*(ϕ-ϕ_s)^2 ] * S(ϕ-ϕ_s,λ)
static void Inflation_V_phi_I_func(arb_t res, const arb_t phi, const slong order, slong prec)
{
    arb_t s,t,w,delta_phi;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(delta_phi);
    
    arb_sub(delta_phi,phi,Inf_Phi_s,prec); //Δϕ=ϕ-ϕ_s
    
    switch ( order )
    {
        case 0: //原函数
            // sqrt_E1*(ϕ-ϕ_s) + η_1*(ϕ-ϕ_s)^2
            arb_mul(s,Inf_sqrt_E1,delta_phi,prec);
            
            arb_sqr(w,delta_phi,prec);
            arb_mul(t,Inf_Eta_1,w,prec);
            
            arb_add(res,s,t,prec);
            break;
        case 1: //一阶导
            // sqrt_E1 + 2*η_1*(ϕ-ϕ_s)
            arb_mul(s,Inf_Eta_1,delta_phi,prec);
            arb_mul_ui(s,s,2,prec);
            
            arb_add(res,s,Inf_sqrt_E1,prec);
            break;
        case 2: //二阶导
            // 2*η_1
            arb_mul_ui(res,Inf_Eta_1,2,prec);
            
            break;
        default:
            printf("src/func/Power_Spectra_numeric/02_various_func.c Inflation_V_phi_I_func 阶数输入有误\n");
            exit(-1);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(delta_phi);
}


//第二段势能，及其导数
//V_0 在总势能那里统一处理
//第二个SR阶段，V_II(ϕ)=[ ΔV + sqrt_E2*(ϕ-ϕ_e) + η_2*(ϕ-ϕ_e)^2 ]*S(-(ϕ-ϕ_e),λ)
static void Inflation_V_phi_II_func(arb_t res, const arb_t phi, const slong order, slong prec)
{
    arb_t s,t,w,delta_phi;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(delta_phi);
    
    arb_sub(delta_phi,phi,Inf_Phi_e,prec); //Δϕ=ϕ-ϕ_e
    
    switch ( order )
    {
        case 0: //原函数
            // ΔV + sqrt_E2*(ϕ-ϕ_e) + η_2*(ϕ-ϕ_e)^2
            arb_mul(s,Inf_sqrt_E2,delta_phi,prec);
            arb_add(s,s,Inf_Delta_V,prec);
            
            arb_sqr(w,delta_phi,prec);
            arb_mul(t,Inf_Eta_2,w,prec);
            
            arb_add(res,s,t,prec);
            
            break;
        case 1: //一阶导
            // sqrt_E2 + 2*η_2*(ϕ-ϕ_e)
            arb_mul(s,Inf_Eta_2,delta_phi,prec);
            arb_mul_ui(s,s,2,prec);
            
            arb_add(res,s,Inf_sqrt_E2,prec);
            
            break;
        case 2: //二阶导
            // 2*η_2
            
            arb_mul_ui(res,Inf_Eta_2,2,prec);
            break;
        default:
            printf("src/func/Power_Spectra_numeric/02_various_func.c Inflation_V_phi_II_func 阶数输入有误\n");
            exit(-1);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(delta_phi);
}


//暴胀子总的势能，及其导数 V(ϕ)
void Inflation_V_phi(arb_t res, const arb_t phi, const slong order, slong prec)
{
    arb_t s,t,w,phi_small,phi_large,s_0_small,s_0_large,s_1_small,s_1_large,s_2_small,s_2_large;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(phi_small);
    arb_init(phi_large);
    arb_init(s_0_small); //原函数用
    arb_init(s_0_large);
    arb_init(s_1_small); //一阶导用
    arb_init(s_1_large);
    arb_init(s_2_small); //二阶导用
    arb_init(s_2_large);
    
    //phi after the step 
    arb_sub(phi_small,phi,Inf_Phi_e,prec); // -(ϕ-ϕ_e)
    arb_neg(phi_small,phi_small); 
    arb_sub(phi_large,phi,Inf_Phi_s,prec); //ϕ-ϕ_s
    
    Inflation_Potential_Smoothing_Function(s_0_small, phi_small, Inf_Lambda, 0, 1, prec); //原函数
    Inflation_Potential_Smoothing_Function(s_0_large, phi_large, Inf_Lambda, 0, 2, prec); //原函数
    
    switch ( order )
    {
        case 0: //原函数
            //V = V0 + V_I(ϕ) * S(ϕ-ϕ_s,λ) + V_II(ϕ) * S(-(ϕ-ϕ_e),λ)
            
            Inflation_V_phi_I_func(t, phi, 0, prec); //原函数
            arb_mul(t,t,s_0_large,prec);
            
            Inflation_V_phi_II_func(s, phi, 0, prec); //原函数
            arb_mul(s,s,s_0_small,prec);
            
            arb_add(s,s,t,prec);
            arb_add(res,s,Inf_V0,prec); //加上最前面的 V_0
            break;
            
        case 1: //一阶导
            // V' = V'_I(ϕ) * S(ϕ-ϕ_s,λ) + V'_II(ϕ) * S(-(ϕ-ϕ_e),λ) 
            //      + V_I(ϕ) * S'(ϕ-ϕ_s,λ)*λ + V_II(ϕ) * S'(-(ϕ-ϕ_e),λ)*(-λ)
            
            Inflation_Potential_Smoothing_Function(s_1_small, phi_small, Inf_Lambda, 1, 1, prec); //一阶导
            Inflation_Potential_Smoothing_Function(s_1_large, phi_large, Inf_Lambda, 1, 2, prec); //一阶导
            
            Inflation_V_phi_I_func(s, phi, 1, prec); //一阶导
            arb_mul(s,s,s_0_large,prec);
            
            Inflation_V_phi_II_func(t, phi, 1, prec); //一阶导
            arb_mul(t,t,s_0_small,prec);
            arb_add(s,s,t,prec);
            
            Inflation_V_phi_I_func(t, phi, 0, prec); //原函数
            arb_mul(t,t,s_1_large,prec);
            arb_mul(t,t,Inf_Lambda,prec);
            arb_add(s,s,t,prec);
            
            Inflation_V_phi_II_func(t, phi, 0, prec); //原函数
            arb_mul(t,t,s_1_small,prec);
            arb_mul(t,t,Inf_Lambda,prec);
            arb_neg(t,t);
            arb_add(res,s,t,prec);
            
            break;
        case 2: //二阶导
            // V'' = V''_I(ϕ) * S(ϕ-ϕ_s,λ) + V''_II(ϕ) * S(-(ϕ-ϕ_e),λ) + V'_I(ϕ) * S'(ϕ-ϕ_s,λ)*λ + V'_II(ϕ) * S'(-(ϕ-ϕ_e),λ)*(-λ)
            //      + V'_I(ϕ) * S'(ϕ-ϕ_s,λ)*λ + V'_II(ϕ) * S'(-(ϕ-ϕ_e),λ)*(-λ) + V_I(ϕ) * S''(ϕ-ϕ_s,λ)*λ^2 + V_II(ϕ) * S''(-(ϕ-ϕ_e),λ)*λ^2
            
            // V'' = V''_I(ϕ) * S(ϕ-ϕ_s,λ) + V''_II(ϕ) * S(-(ϕ-ϕ_e),λ)
            //       + 2*V'_I(ϕ) * S'(ϕ-ϕ_s,λ)*λ + 2*V'_II(ϕ) * S'(-(ϕ-ϕ_e),λ)*(-λ)
            //       + V_I(ϕ) * S''(ϕ-ϕ_s,λ)*λ^2 + V_II(ϕ) * S''(-(ϕ-ϕ_e),λ)*λ^2
            
            
            Inflation_Potential_Smoothing_Function(s_1_small, phi_small, Inf_Lambda, 1, 1, prec); //一阶导
            Inflation_Potential_Smoothing_Function(s_1_large, phi_large, Inf_Lambda, 1, 2, prec); //一阶导
            
            Inflation_Potential_Smoothing_Function(s_2_small, phi_small, Inf_Lambda, 2, 1, prec); //二阶导
            Inflation_Potential_Smoothing_Function(s_2_large, phi_large, Inf_Lambda, 2, 2, prec); //二阶导
            
            Inflation_V_phi_I_func(s, phi, 2, prec); //二阶导
            arb_mul(s,s,s_0_large,prec);
            
            Inflation_V_phi_II_func(t, phi, 2, prec); //二阶导
            arb_mul(t,t,s_0_small,prec);
            arb_add(s,s,t,prec);
            
            Inflation_V_phi_I_func(t, phi, 1, prec); //一阶导
            arb_mul(t,t,s_1_large,prec);
            arb_mul_ui(t,t,2,prec);
            arb_mul(t,t,Inf_Lambda,prec);
            arb_add(s,s,t,prec);
            
            Inflation_V_phi_II_func(t, phi, 1, prec); //一阶导
            arb_mul(t,t,s_1_small,prec);
            arb_mul_ui(t,t,2,prec);
            arb_mul(t,t,Inf_Lambda,prec);
            arb_neg(t,t);
            arb_add(s,s,t,prec);
            
            Inflation_V_phi_I_func(t, phi, 0, prec); //原函数
            arb_mul(t,t,s_2_large,prec);
            arb_mul(t,t,Inf_Lambda,prec);
            arb_mul(t,t,Inf_Lambda,prec);
            arb_add(s,s,t,prec);
            
            Inflation_V_phi_II_func(t, phi, 0, prec); //原函数
            arb_mul(t,t,s_2_small,prec);
            arb_mul(t,t,Inf_Lambda,prec);
            arb_mul(t,t,Inf_Lambda,prec);
            arb_add(res,s,t,prec);
            
            break;
        default:
            printf("src/func/Power_Spectra_numeric/02_various_func.c Inflation_V_phi 阶数输入有误\n");
            exit(-1);
    }
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(phi_small);
    arb_clear(phi_large);
    arb_clear(s_0_small);
    arb_clear(s_0_large);
    arb_clear(s_1_small);
    arb_clear(s_1_large);
    arb_clear(s_2_small);
    arb_clear(s_2_large);
}

