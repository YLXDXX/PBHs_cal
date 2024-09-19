#include "Inflation_odes_slove.h"
#include <stdlib.h> 

//暴胀模型各种初始参数设定：势能参数、慢滚参数等

// model parameters
arb_t Inf_Mpl,Inf_Mpl_2,Inf_Mpl_4,Inf_Mpl_6; //Planck质量极其幂次
arb_t Inf_Delta_V,Inf_V0,Inf_Phi_e;
arb_t Inf_Epsilon_1,Inf_sqrt_E1,Inf_Epsilon_2,Inf_sqrt_E2;
arb_t Inf_Eta_1,Inf_Eta_2;
arb_t Inf_Lambda,Inf_Delta_Phi_usr,Inf_Phi_s;
static arb_ptr Inf_ODEs_model_parameters_help; //仅用于参数初始化


// Set model parameters
void Inflation_set_model_parameters(slong prec)
{
    if(Inf_ODEs_model_parameters_help!=NULL)
    {
        return;
    }
    
    Inf_ODEs_model_parameters_help=_arb_vec_init(1);
    
    arb_init(Inf_Mpl);
    arb_init(Inf_Mpl_2);
    arb_init(Inf_Mpl_4);
    arb_init(Inf_Mpl_6);
    arb_init(Inf_Delta_V);
    arb_init(Inf_V0);
    arb_init(Inf_Phi_e);
    
    arb_init(Inf_Epsilon_1);
    arb_init(Inf_sqrt_E1);
    arb_init(Inf_Epsilon_2);
    arb_init(Inf_sqrt_E2);
    arb_init(Inf_Eta_1);
    arb_init(Inf_Eta_2);
    arb_init(Inf_Lambda);
    arb_init(Inf_Delta_Phi_usr);
    arb_init(Inf_Phi_s);
    
    // Planck mass 
    arb_set_str(Inf_Mpl,"1",prec); //Planck 质量 Mpl
    arb_pow_ui(Inf_Mpl_2,Inf_Mpl,2,prec); //Mpl^2
    arb_pow_ui(Inf_Mpl_4,Inf_Mpl,4,prec); //Mpl^4
    arb_pow_ui(Inf_Mpl_6,Inf_Mpl,6,prec); //Mpl^6
    
    // Height of upward step (1.264*10^-15*Mpl^4 in previous version)
    arb_set_str(Inf_Delta_V,"7.6342e-13",prec); //势能 ΔV = 7.6342e-13 * Mpl**4
    //arb_set_str(Inf_Delta_V,"5.97388e-13",prec); //势能 ΔV = 5.97388e-13 * Mpl**4
    //arb_set_str(Inf_Delta_V,"5.94104e-13",prec); //势能 ΔV = 5.94104e-13 * Mpl**4
    arb_mul(Inf_Delta_V,Inf_Delta_V,Inf_Mpl_4,prec);
    
    // Potential at initial stage (previous version: 10^-10)
    arb_set_str(Inf_V0,"7e-10",prec); //势能 V_0 = 7e-10 * Mpl**4
    arb_mul(Inf_V0,Inf_V0,Inf_Mpl_4,prec);
    
    // Field value at step
    arb_mul_ui(Inf_Phi_e,Inf_Mpl,5,prec); // ϕ_e = 5 * Mpl 
    
    // First slow-roll parameter in the first slow-roll stage
    arb_set_str(Inf_Epsilon_1,"1.25e-21",prec); // ε_1 = 1.25e-21 * Mpl**6 
    arb_mul(Inf_Epsilon_1,Inf_Epsilon_1,Inf_Mpl_6,prec);
    
    arb_mul_ui(Inf_sqrt_E1,Inf_Epsilon_1,2,prec); // sqrt(2*ε_1) = sqrt(2.*Epsilon_1)
    arb_sqrt(Inf_sqrt_E1,Inf_sqrt_E1,prec); 
    
    
    // First slow-roll parameter in the second slow-roll stage
    arb_set_str(Inf_Epsilon_2,"1.25e-21",prec); // ε_2 = 1.8e-25 * Mpl**6
    arb_mul(Inf_Epsilon_2,Inf_Epsilon_2,Inf_Mpl_6,prec);
    
    arb_mul_ui(Inf_sqrt_E2,Inf_Epsilon_2,2,prec); // sqrt(2*ε_2) = sqrt(2.*Epsilon_2)
    arb_sqrt(Inf_sqrt_E2,Inf_sqrt_E2,prec); 
    
    
    // Second slow-roll parameter in the first slow-roll stage
    arb_set_str(Inf_Eta_1,"-1e-11",prec); // η_1 = -1e-11 * Mpl**2
    arb_mul(Inf_Eta_1,Inf_Eta_1,Inf_Mpl_2,prec);
    
    // Second slow-roll parameter in the second slow-roll stage
    arb_set_str(Inf_Eta_2,"0",prec); // η_2 = 0 * Mpl**2
    arb_mul(Inf_Eta_2,Inf_Eta_2,Inf_Mpl_2,prec);
    
    
    // Used to adjust the steepness of the Tanh function
    arb_set_str(Inf_Lambda,"1e2",prec); // λ = 1e2 / Mpl
    //arb_set_str(Inf_Lambda,"5e3",prec); // λ = 5e3 / Mpl
    //arb_set_str(Inf_Lambda,"5e4",prec); // λ = 5e4 / Mpl
    arb_div(Inf_Lambda,Inf_Lambda,Inf_Mpl,prec);
    
    
    // Field range of ultra slow roll
    arb_set_str(Inf_Delta_Phi_usr,"0",prec); // Δϕ = 0.02 * Mpl
    arb_mul(Inf_Delta_Phi_usr,Inf_Delta_Phi_usr,Inf_Mpl,prec);
    
    // Field value when ultra slow roll starts
    arb_add(Inf_Phi_s,Inf_Phi_e,Inf_Delta_Phi_usr,prec); // ϕ_s = ϕ_e + Δϕ （ Phi_s = Phi_e + Delta_Phi_usr ）
}


