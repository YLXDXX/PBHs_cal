#include "Inflation_odes_slove.h"
#include <stdlib.h>
#include <string.h> 


Interp_coe_t G_Interp_coe_0;
Interp_coe_t G_Interp_coe_1;
Interp_coe_t G_Interp_coe_2;
Interp_coe_t G_Interp_coe_3;
Interp_coe_t G_Interp_coe_4;
arb_t G_fourier_k;

//扰动求解


// model parameters
static arb_t Mpl,Mpl_2,Mpl_4,Mpl_6,Delta_V,V0,Phi_e;
static arb_t Epsilon_1,sqrt_E1,Epsilon_2,sqrt_E2,Eta_1,Eta_2,Lambda,Delta_Phi_usr,Phi_s;
static arb_ptr ttttttttt;

// Set model parameters
void Set_ODEs_model_parameters(slong prec)
{
    if(ttttttttt!=NULL)
    {
        return;
    }
    
    ttttttttt=_arb_vec_init(1);
    
    arb_init(Mpl);
    arb_init(Mpl_2);
    arb_init(Mpl_4);
    arb_init(Mpl_6);
    arb_init(Delta_V);
    arb_init(V0);
    arb_init(Phi_e);
    
    arb_init(Epsilon_1);
    arb_init(sqrt_E1);
    arb_init(Epsilon_2);
    arb_init(sqrt_E2);
    arb_init(Eta_1);
    arb_init(Eta_2);
    arb_init(Lambda);
    arb_init(Delta_Phi_usr);
    arb_init(Phi_s);
    
    // Planck mass 
    arb_set_str(Mpl,"1",prec); 
    arb_pow_ui(Mpl_2,Mpl,2,prec); //Mpl^2
    arb_pow_ui(Mpl_4,Mpl,4,prec); //Mpl^4
    arb_pow_ui(Mpl_6,Mpl,6,prec); //Mpl^6
    
    // Height of upward step (1.264*10^-15*Mpl^4 in previous version)
    arb_set_str(Delta_V,"6.3e-16",prec); //6.3e-16 * Mpl**4
    arb_mul(Delta_V,Delta_V,Mpl_4,prec);
    
    // Potential at initial stage (previous version: 10^-10)
    arb_set_str(V0,"7e-10",prec); //7e-10 * Mpl**4
    arb_mul(V0,V0,Mpl_4,prec);
    
    // Field value at step
    arb_mul_ui(Phi_e,Mpl,5,prec); // 5 * Mpl 
    
    // First slow-roll parameter in the first slow-roll stage
    arb_set_str(Epsilon_1,"1.25e-21",prec); //1.25e-21 * Mpl**6 
    arb_mul(Epsilon_1,Epsilon_1,Mpl_6,prec);
    
    arb_mul_ui(sqrt_E1,Epsilon_1,2,prec); //sqrt(2.*Epsilon_1)
    arb_sqrt(sqrt_E1,sqrt_E1,prec); 
    
    
    // First slow-roll parameter in the second slow-roll stage
    arb_set_str(Epsilon_2,"1.8e-25",prec); //1.8e-25 * Mpl**6
    arb_mul(Epsilon_2,Epsilon_2,Mpl_6,prec);
    
    arb_mul_ui(sqrt_E2,Epsilon_2,2,prec); //sqrt(2.*Epsilon_2)
    arb_sqrt(sqrt_E2,sqrt_E2,prec); 
    
    
    // Second slow-roll parameter in the first slow-roll stage
    arb_set_str(Eta_1,"-1e-11",prec); //-1e-11 * Mpl**2
    arb_mul(Eta_1,Eta_1,Mpl_2,prec);
    
    // Second slow-roll parameter in the second slow-roll stage
    arb_set_str(Eta_2,"0",prec); //0 * Mpl**2
    arb_mul(Eta_2,Eta_2,Mpl_2,prec);
    
    
    // Used to adjust the steepness of the Tanh function
    arb_set_str(Lambda,"5e4",prec); //5e4 / Mpl
    arb_div(Lambda,Lambda,Mpl,prec);
    
    
    // Field range of ultra slow roll
    arb_set_str(Delta_Phi_usr,"0.02",prec); //0.02 * Mpl
    arb_mul(Delta_Phi_usr,Delta_Phi_usr,Mpl,prec);
    
    // Field value when ultra slow roll starts
    arb_add(Phi_s,Phi_e,Delta_Phi_usr,prec); //Phi_e + Delta_Phi_usr
}


//求解背景扰动方程所需函数

void Func_Smoothing_Step_Function(arb_t res, const arb_t x, const arb_t smoothness, slong prec)
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
        //(t**2) * (3 - 2.*t)
        arb_sqr(s,t,prec);
        arb_mul_ui(w,t,2,prec);
        arb_neg(w,w);
        arb_add_ui(w,w,3,prec);
        arb_mul(res,s,w,prec);
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

void Func_S_prime(arb_t res, const arb_t x, const arb_t smoothness, slong prec)
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
        //6.*(t) * (1. - t)
        arb_mul_ui(s,t,6,prec);
        arb_neg(w,t);
        arb_add_ui(w,w,1,prec);
        arb_mul(res,s,w,prec);
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


void Func_S_pp(arb_t res, const arb_t x, const arb_t smoothness, slong prec)
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


void Func_V_phi(arb_t res, const arb_t phi, slong prec)
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
    arb_sub(phi_small,phi,Phi_e,prec);
    arb_sub(phi_large,phi,Phi_s,prec);
    
    arb_neg(w,phi_small);
    Func_Smoothing_Step_Function(smooth_small, w, Lambda, prec);
    Func_Smoothing_Step_Function(smooth_large, phi_large, Lambda, prec);
    
    //V = (V0 + (Delta_V + sqrt_E2 * phi_small + Eta_2 * phi_small**2) * smooth_small + 
    //     sqrt_E1 * phi_large * smooth_large)
    
    //中间部分
    //(Delta_V + sqrt_E2 * phi_small + Eta_2 * phi_small**2) * smooth_small
    arb_mul(s,sqrt_E2,phi_small,prec);
    arb_sqr(t,phi_small,prec);
    arb_mul(t,t,Eta_2,prec);
    arb_add(s,s,t,prec);
    arb_add(s,s,Delta_V,prec);
    arb_mul(s,s,smooth_small,prec);
    
    arb_add(s,s,V0,prec); //加上最前面的 V0
    
    //未尾部分
    //sqrt_E1 * phi_large * smooth_large
    
    arb_mul(t,sqrt_E1,phi_large,prec);
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
void Func_V_phi_p(arb_t res, const arb_t phi, slong prec)
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
    arb_sub(phi_small,phi,Phi_e,prec);
    arb_sub(phi_large,phi,Phi_s,prec);
    
    arb_neg(w,phi_small);
    Func_Smoothing_Step_Function(smooth_small, w, Lambda, prec);
    Func_Smoothing_Step_Function(smooth_large, phi_large, Lambda, prec);
    
    //((Delta_V + sqrt_E2 * phi_small + Eta_2 * phi_small**2) * (-Lambda * S_prime(-phi_small, Lambda)) +
    //      (sqrt_E2 + Eta_2 * phi_small * 2.) * smooth_small + 
    //      sqrt_E1 * smooth_large + sqrt_E1 * phi_large * S_prime(phi_large, Lambda) * Lambda)
    
    //前面
    //(Delta_V + sqrt_E2 * phi_small + Eta_2 * phi_small**2) * (-Lambda * S_prime(-phi_small, Lambda))
    arb_mul(s,sqrt_E2,phi_small,prec);
    arb_sqr(t,phi_small,prec);
    arb_mul(t,t,Eta_2,prec);
    arb_add(s,s,t,prec);
    arb_add(s,s,Delta_V,prec);
    
    arb_neg(t,phi_small);
    Func_S_prime(w, t, Lambda, prec);
    arb_mul(w,w,Lambda,prec);
    arb_neg(w,w);
    arb_mul(s,s,w,prec);
    
    //中间
    //(sqrt_E2 + Eta_2 * phi_small * 2.) * smooth_small + sqrt_E1 * smooth_large
    arb_mul(t,Eta_2,phi_small,prec);
    arb_mul_ui(t,t,2,prec);
    arb_add(t,t,sqrt_E2,prec);
    arb_mul(t,t,smooth_small,prec);
    arb_add(s,s,t,prec);
    
    arb_mul(t,sqrt_E1,smooth_large,prec);
    arb_add(s,s,t,prec);
    
    //最后
    //sqrt_E1 * phi_large * S_prime(phi_large, Lambda) * Lambda
    
    arb_mul(t,sqrt_E1,phi_large,prec);
    Func_S_prime(w, phi_large, Lambda, prec);
    arb_mul(t,t,w,prec);
    arb_mul(t,t,Lambda,prec);
    
    arb_add(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(phi_small);
    arb_clear(phi_large);
    arb_clear(smooth_small);
    arb_clear(smooth_large);
}

void Func_V_phi_pp(arb_t res, const arb_t phi, slong prec)
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
    arb_sub(phi_small,phi,Phi_e,prec);
    arb_sub(phi_large,phi,Phi_s,prec);
    
    arb_neg(w,phi_small);
    Func_Smoothing_Step_Function(smooth_small, w, Lambda, prec);
    Func_Smoothing_Step_Function(smooth_large, phi_large, Lambda, prec);
    
    // ((Delta_V + sqrt_E2 * phi_small + Eta_2 * phi_small**2) * ((Lambda**2) * S_pp(-phi_small, Lambda)) +
    //      (sqrt_E2 + Eta_2 * phi_small*2.) * (-Lambda * S_prime(-phi_small, Lambda)) +
    //      ( Eta_2 * 2.) * smooth_small + 
    //      (sqrt_E2 + Eta_2 * phi_small * 2.) * S_prime(-phi_small, Lambda)*(-Lambda ) +
    //      sqrt_E1 * S_prime(phi_large, Lambda)*(Lambda) + 
    //      sqrt_E1 * S_prime(phi_large, Lambda) * Lambda + 
    //      sqrt_E1 * phi_large * S_pp(phi_large, Lambda) * (Lambda**2.))
          
    //前半部分
    //(Delta_V + sqrt_E2 * phi_small + Eta_2 * phi_small**2) * ((Lambda**2) * S_pp(-phi_small, Lambda))
    arb_mul(s,sqrt_E2,phi_small,prec);
    arb_sqr(t,phi_small,prec);
    arb_mul(t,t,Eta_2,prec);
    arb_add(s,s,t,prec);
    arb_add(s,s,Delta_V,prec);
    
    arb_neg(t,phi_small);
    Func_S_pp(w, t, Lambda, prec);
    arb_sqr(t,Lambda,prec);
    arb_mul(t,t,w,prec);
    arb_mul(s,s,t,prec);
    
    //中间 (sqrt_E2 + Eta_2 * phi_small*2.) * (-Lambda * S_prime(-phi_small, Lambda))
    arb_neg(t,phi_small);
    Func_S_prime(w, t, Lambda, prec);
    arb_mul(w,w,Lambda,prec);
    arb_neg(w,w);
    
    arb_mul(t,Eta_2,phi_small,prec);
    arb_mul_si(t,t,2,prec);
    arb_add(t,t,sqrt_E2,prec);
    
    arb_mul(w,t,w,prec);
    arb_add(s,s,w,prec);
    
    //中间 ( Eta_2 * 2.) * smooth_small
    arb_mul_si(t,Eta_2,2,prec);
    arb_mul(t,t,smooth_small,prec);
    arb_add(s,s,t,prec);
    
    //中间 (sqrt_E2 + Eta_2 * phi_small * 2.) * S_prime(-phi_small, Lambda)*(-Lambda )
    //这里，与前面一样
    arb_add(s,s,w,prec);
    
    //中间 sqrt_E1 * S_prime(phi_large, Lambda)*(Lambda)
    Func_S_prime(w, phi_large, Lambda, prec);
    arb_mul(w,w,sqrt_E1,prec);
    arb_mul(w,w,Lambda,prec);
    
    arb_add(s,s,w,prec);
    
    //中间 sqrt_E1 * S_prime(phi_large, Lambda) * Lambda
    //这里，与前面一样
    arb_add(s,s,w,prec);
    
    //最后 sqrt_E1 * phi_large * S_pp(phi_large, Lambda) * (Lambda**2.)
    
    arb_mul(t,sqrt_E1,phi_large,prec);
    Func_S_pp(w, phi_large, Lambda, prec);
    arb_mul(t,t,w,prec);
    arb_sqr(w,Lambda,prec);
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



// define the ODE system
int Func_background_coupled_odes(arb_ptr yp, const arb_t t, const arb_ptr y, const slong dim,
                           void* param, const slong order, slong prec)
{
    Set_ODEs_model_parameters(prec);
    
    arb_t s,w;
    arb_t phi_dot,phi,tau,N,H,V,dV_dphi;
    
    arb_init(s);
    arb_init(w);
    arb_init(phi_dot);
    arb_init(phi);
    arb_init(tau);
    arb_init(N);
    arb_init(H);
    arb_init(V);
    arb_init(dV_dphi);
    
    // y = {phi_dot, phi, tau, N} 
    // yp= {dphi_dot_dt, dphi_dt, dtau_dt, dN_dt}
    
    arb_set(phi_dot,y);
    arb_set(phi,y+1);
    arb_set(tau,y+2);
    arb_set(N,y+3);
    
    // Calculate the potential and its derivative
    Func_V_phi(V,phi,prec); // V = V_phi(phi)
    Func_V_phi_p(dV_dphi,phi,prec); // dV_dphi = V_phi_p(phi)
    
    // Calculate H based on the current values of phi and phi_dot
    
    //H = np.sqrt((1./6.) * phi_dot**2 + V / 3.)
    arb_sqr(s,phi_dot,prec);
    arb_div_ui(s,s,6,prec);
    arb_div_ui(w,V,3,prec);
    arb_add(s,s,w,prec);
    arb_sqrt(H,s,prec);
    
    // ODEs
    //dphi_dot_dt = - (3. * H * phi_dot + dV_dphi)
    //dphi_dt = phi_dot
    //dtau_dt = 1. / np.exp(N)
    //dN_dt = H
    
    //第一个
    arb_mul(s,H,phi_dot,prec);
    arb_mul_si(s,s,3,prec);
    arb_add(s,s,dV_dphi,prec);
    arb_neg(yp,s);
    
    //第二个
    arb_set(yp+1,phi_dot);
    
    //第三个
    arb_exp(s,N,prec);
    arb_inv(yp+2,s,prec);
    
    //第四个
    arb_set(yp+3,H);
    
    arb_clear(s);
    arb_clear(w);
    arb_clear(phi_dot);
    arb_clear(phi);
    arb_clear(tau);
    arb_clear(N);
    arb_clear(H);
    arb_clear(V);
    arb_clear(dV_dphi);
    
    return 0;
}



void Func_m_eff(arb_t res, const arb_t phi_dot, const arb_t phi_ddot,
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
    
    arb_inv(t,Mpl,prec);
    arb_sqr(t,t,prec);
    arb_mul(s,s,t,prec);
    
    arb_sub(res,V_phi_pp,s,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
}

int Func_perturbation_phi_odes(arb_ptr yp, const arb_t t, const arb_ptr y, const slong dim,
                                 void* param, const slong order, slong prec)
{
    //注意到，这里扰动的方程需要用到背景的解
    //这里将背景制的微分方程多点解，通过参数传入
    ODEs_DOPRI54_dense_t d_out;
    d_out=param;
    
    arb_t s,w,q;
    arb_t Qphi_Real,Qphi_Real_dot,Qphi_Imag,Qphi_Imag_dot;
    arb_t phi,phi_dot,V,dV_dphi,d2V_dphi2;
    arb_t H,H_dot,N,phi_ddot,fk;
    
    arb_init(s);
    arb_init(w);
    arb_init(q);
    arb_init(Qphi_Real);
    arb_init(Qphi_Real_dot);
    arb_init(Qphi_Imag);
    arb_init(Qphi_Imag_dot);
    
    arb_init(phi);
    arb_init(phi_dot);
    arb_init(V);
    arb_init(dV_dphi);
    arb_init(d2V_dphi2);
    
    arb_init(H);
    arb_init(H_dot);
    arb_init(N);
    arb_init(phi_ddot);
    arb_init(fk);
    
    //fk define the fourier k modes
    arb_set(fk,G_fourier_k);
    
    //{Qphi_Real, Qphi_Real_dot, Qphi_Imag, Qphi_Imag_dot} = y
    arb_set(Qphi_Real,y);
    arb_set(Qphi_Real_dot,y+1);
    arb_set(Qphi_Imag,y+2);
    arb_set(Qphi_Imag_dot,y+3);
    
    
    // Calculate the potential and its derivative
    //背景解 phi = phi_interp(t)
    Interpolation_fit_func_odes_DOPRI54(phi, t, d_out, 1, prec);
    
    //背景解 phi_dot = phi_dot_interp(t)
    Interpolation_fit_func_odes_DOPRI54(phi_dot, t, d_out, 0, prec);
    
    Func_V_phi(V,phi,prec); //V = V_phi(phi)
    Func_V_phi_p(dV_dphi,phi,prec); //dV_dphi = V_phi_p(phi)
    Func_V_phi_pp(d2V_dphi2,phi,prec); //d2V_dphi2 = V_phi_pp(phi) //The second derivative of V(\phi)
    
    //背景解 H = H_interp(t), 这里不用插值
    //H = np.sqrt((1./6.) * phi_dot**2 + V / 3.)
    arb_sqr(s,phi_dot,prec);
    arb_div_ui(s,s,6,prec);
    arb_div_ui(w,V,3,prec);
    arb_add(s,s,w,prec);
    arb_sqrt(H,s,prec);
    
    //H_dot= -0.5*(phi_dot**2.)
    arb_sqr(w,phi_dot,prec);
    arb_div_si(H_dot,w,-2,prec);
    
    //背景解 N = N_interp(t)
    Interpolation_fit_func_odes_DOPRI54(N, t, d_out, 3,  prec);
    
    //phi_ddot = - (3. * H * phi_dot + dV_dphi)
    arb_mul(w,H,phi_dot,prec);
    arb_mul_si(w,w,3,prec);
    arb_add(w,w,dV_dphi,prec);
    arb_neg(phi_ddot,w);
    
    
    // ODEs
    // yp = { dQphi_Real_dt, dQphi_Real_dot_dt, dQphi_Imag_dt, dQphi_Imag_dot_dt }
    
    //第一个方程
    //dQphi_Real_dt = Qphi_Real_dot
    arb_set(yp,Qphi_Real_dot);
    
    //第二个方程
    //dQphi_Real_dot_dt = -(3*H*Qphi_Real_dot + (fk**2/np.exp(2.*N))*Qphi_Real + m_eff(phi_dot, phi_ddot, H, H_dot, d2V_dphi2)*Qphi_Real)
    //3*H*Qphi_Real_dot
    arb_mul(s,H,Qphi_Real_dot,prec);
    arb_mul_ui(s,s,3,prec);
    
    //(fk**2/np.exp(2.*N))*Qphi_Real
    arb_sqr(w,fk,prec);
    arb_mul_ui(q,N,2,prec);
    arb_exp(q,q,prec);
    arb_div(w,w,q,prec);
    arb_mul(w,w,Qphi_Real,prec);
    arb_add(s,s,w,prec);
    
    //m_eff(phi_dot, phi_ddot, H, H_dot, d2V_dphi2)*Qphi_Real
    Func_m_eff(w,phi_dot, phi_ddot, H, H_dot, d2V_dphi2,prec);
    arb_mul(w,w,Qphi_Real,prec);
    arb_add(s,s,w,prec);
    
    arb_neg(yp+1,s);
    
    //第三个方程
    //dQphi_Imag_dt = Qphi_Imag_dot
    arb_set(yp+2,Qphi_Imag_dot);
    
    //第四个方程
    //dQphi_Imag_dot_dt = -(3*H*Qphi_Imag_dot + (fk**2/np.exp(2.*N))*Qphi_Imag + m_eff(phi_dot, phi_ddot, H, H_dot, d2V_dphi2)*Qphi_Imag)
    
    //3*H*Qphi_Imag_dot
    arb_mul(s,H,Qphi_Imag_dot,prec);
    arb_mul_ui(s,s,3,prec);
    
    //(fk**2/np.exp(2.*N))*Qphi_Imag
    arb_sqr(w,fk,prec);
    arb_mul_ui(q,N,2,prec);
    arb_exp(q,q,prec);
    arb_div(w,w,q,prec);
    arb_mul(w,w,Qphi_Imag,prec);
    arb_add(s,s,w,prec);
    
    //m_eff(phi_dot, phi_ddot, H, H_dot, d2V_dphi2)*Qphi_Imag
    Func_m_eff(w,phi_dot, phi_ddot, H, H_dot, d2V_dphi2,prec);
    arb_mul(w,w,Qphi_Imag,prec);
    arb_add(s,s,w,prec);
    
    arb_neg(yp+3,s);
    
    
    arb_clear(s);
    arb_clear(w);
    arb_clear(q);
    arb_clear(Qphi_Real);
    arb_clear(Qphi_Real_dot);
    arb_clear(Qphi_Imag);
    arb_clear(Qphi_Imag_dot);
    
    arb_clear(phi);
    arb_clear(phi_dot);
    arb_clear(V);
    arb_clear(dV_dphi);
    arb_clear(d2V_dphi2);
    
    arb_clear(H);
    arb_clear(H_dot);
    arb_clear(N);
    arb_clear(phi_ddot);
    arb_clear(fk);
    
    return 0;
}


struct Get_time_k_enter_structure {
    ODEs_DOPRI54_dense_t d_out;
    arb_t ln_a_i;
    arb_t ln_k;
};


int interior_Func_get_time_k_enter(arb_t res, const arb_t ln_t,
                                   void * pass, const slong order,  slong prec)
{
    arb_t s,t,w,N,ln_H,phi,phi_dot,V;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(N);
    arb_init(ln_H);
    arb_init(phi);
    arb_init(phi_dot);
    arb_init(V);
    
    struct Get_time_k_enter_structure* param;
    param=pass;
    
    
    //进入视界条件 k=aH ⇒ ln(k)=N+ln(a_i)+ln(H)
    arb_exp(t,ln_t,prec); //恢复线性时间
    
    //背景解 N = N_interp(t)
    Interpolation_fit_func_odes_DOPRI54(N, t, param->d_out, 3,  prec);
    
    //背景解 phi = phi_interp(t)
    Interpolation_fit_func_odes_DOPRI54(phi, t, param->d_out, 1, prec);
    
    //背景解 phi_dot = phi_dot_interp(t)
    Interpolation_fit_func_odes_DOPRI54(phi_dot, t, param->d_out, 0, prec);
    
    Func_V_phi(V,phi,prec); //V = V_phi(phi)
    
    //背景解 H = H_interp(t), 这里不用插值
    //H = np.sqrt((1./6.) * phi_dot**2 + V / 3.)
    arb_sqr(s,phi_dot,prec);
    arb_div_ui(s,s,6,prec);
    arb_div_ui(w,V,3,prec);
    arb_add(s,s,w,prec);
    arb_sqrt(ln_H,s,prec); //得到 H
    arb_log(ln_H,ln_H,prec); //取对数
    
    //ln(k)=N+ln(a_i)+ln(H) ⇒ N+ln(a_i)+ln(H)-ln(k)=0
    
    arb_add(s,N,param->ln_a_i,prec);
    arb_add(s,s,ln_H,prec);
    
    //arb_set_str(w,"1.2",prec); //添加系数，将时间适当缩小
    //arb_mul(s,s,w,prec);
    
    arb_sub(res,s,param->ln_k,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(N);
    arb_clear(ln_H);
    arb_clear(phi);
    arb_clear(phi_dot);
    arb_clear(V);
    return 0;
}

//得到暴胀中，模式 k 进入视界的时间及此时的 N
//这里的 k 没取对数，以实际 k 输入
//这里，输出的时间 t 没取对数值
void Func_get_time_k_enter(arb_t t, arb_t N, const arb_t k, const arb_t t_a, const arb_t t_b, //在区间 [t_a, t_b] 内找根
                           const arb_t a_i, const ODEs_DOPRI54_dense_t d_out, // a_i 为尺度因子的初始值
                           slong prec)
{
    //进入视界条件 k=aH ⇒ ln(k)=N+ln(a_i)+ln(H)
    arb_t s,ln_ta,ln_tb,error;
    
    arb_init(s);
    arb_init(ln_ta);
    arb_init(ln_tb);
    arb_init(error);
    
    
    struct Get_time_k_enter_structure* param; //用来传递参数
    
    param= (struct Get_time_k_enter_structure*)calloc(1,sizeof(struct Get_time_k_enter_structure)); //分配内存
    
    arb_init(param->ln_a_i); //使用前初始化其中变量
    arb_init(param->ln_k);
    param->d_out=d_out;
    
    arb_log(param->ln_a_i,a_i,prec);
    arb_log(param->ln_k,k,prec);
    
    //这里时间 t 会非常大，为了求根方便，我们用对数值计算
    if( arb_is_zero(t_a) ) //区间端点可能会有零的情况
    {
        arb_set_str(s,"1E-50",prec);
        arb_log(ln_ta,s,prec);
    }else
    {
        arb_log(ln_ta,t_a,prec);
    }
    
    if( arb_is_zero(t_b) )
    {
        arb_set_str(s,"1E-50",prec);
        arb_log(ln_tb,s,prec);
    }else
    {
        arb_log(ln_tb,t_b,prec);
    }
    
    
    
    arb_set_ui(error,1E5); //1E-5
    arb_inv(error,error,prec);
    
    //需要利用求根方法来得到 k 进入视界的时间
    Find_interval_root(s, interior_Func_get_time_k_enter, param, 0,
                       ln_ta, ln_tb, error,
                       10,Root_C_Max,
                       prec);
    arb_exp(t,s,prec);
    
    //求出此时的 N，可以通过此来计算 当前 a， ln(a)=N+ln(a_i) 
    Interpolation_fit_func_odes_DOPRI54(s, t, param->d_out, 3,  prec);
    arb_set(N,s);
    
    arb_clear(s);
    arb_clear(ln_ta);
    arb_clear(ln_tb);
    arb_clear(error);
    
    arb_clear(param->ln_a_i); //释放结构体内存
    arb_clear(param->ln_k);
    free(param);
}
