#include "Inflation_odes_slove.h"
#include <stdlib.h>
#include <string.h> 
#include <omp.h>


//背景的微分方程
// define the ODE system
int Inflation_background_coupled_odes(arb_ptr yp, const arb_t t, const arb_ptr y, const slong dim,
                           void* param, const slong order, slong prec)
{
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
    Inflation_V_phi(V,phi,0,prec); //原函数 V = V_phi(phi)
    Inflation_V_phi(dV_dphi,phi,1,prec); //一阶导 dV_dphi = V_phi_p(phi)
    
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



//扰动微分方程
int Inflation_perturbation_phi_odes(arb_ptr yp, const arb_t t, const arb_ptr y, const slong dim,
                                 void* param, const slong order, slong prec)
{
    //注意到，这里扰动的方程需要用到背景的解
    //这里将背景制的微分方程多点解，通过参数传入
    //计算所需的 Fourier 模式也由参数传入
    //传参需使用结构体 Inflation_perturbation_ODEs_param_pass_t
    
    Inflation_perturb_ODEs_param_t param_p; //参数传递
    param_p=param;
    
    Inflation_dense_t d_out; //dense ouput 结果
    d_out=param_p->d_out;
    
    
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
    arb_set(fk,param_p->fourier_k); //由参数传入
    
    //{Qphi_Real, Qphi_Real_dot, Qphi_Imag, Qphi_Imag_dot} = y
    arb_set(Qphi_Real,y);
    arb_set(Qphi_Real_dot,y+1);
    arb_set(Qphi_Imag,y+2);
    arb_set(Qphi_Imag_dot,y+3);
    
    
    // Calculate the potential and its derivative
    //背景解 phi = phi_interp(t)
    Inflation_interp_fit_func_odes(phi, t, d_out, 1, prec);
    
    //背景解 phi_dot = phi_dot_interp(t)
    Inflation_interp_fit_func_odes(phi_dot, t, d_out, 0, prec);
    
    Inflation_V_phi(V,phi,0,prec); //原函数 V = V_phi(phi)
    Inflation_V_phi(dV_dphi,phi,1,prec); //一阶导 dV_dphi = V_phi_p(phi)
    Inflation_V_phi(d2V_dphi2,phi,2,prec); //二阶导 d2V_dphi2 = V_phi_pp(phi) //The second derivative of V(\phi)
    
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
    Inflation_interp_fit_func_odes(N, t, d_out, 3,  prec);
    
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
    Inflation_m_eff(w,phi_dot, phi_ddot, H, H_dot, d2V_dphi2,prec);
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
    Inflation_m_eff(w,phi_dot, phi_ddot, H, H_dot, d2V_dphi2,prec);
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











