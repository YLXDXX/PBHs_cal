#ifndef __PBHS_GW_INDUCED_THIRD_POWER_SPECTRA__   /* Include guard */
#define __PBHS_GW_INDUCED_THIRD_POWER_SPECTRA__

#include "../../new_type.h"
#include "../../phy_constant.h"
#include "../../math_method.h"
#include "../../general.h"
#include "../../other.h"

//该结构体用于与 cuba 库中积分的交互
struct User_Data
{
    arb_ptr interval_a;
    arb_ptr interval_b;
    
    func_NULL *func; //函数指针数组的指针
    
    slong prec;
    void* param;
    slong order;
};


//方法三 2305.19950 中高斯情况的辅助函数
void J_2_u_1_v_1_limited(arb_t res, const arb_t u_1, const arb_t v_1, slong prec);
void J_u_1_v_1_and_u_2_v_2_limited(arb_t res, const arb_t u_1, const arb_t v_1, const arb_t u_2, const arb_t v_2, slong prec);


//多维积分内部函数
int interior_GW_current_energy_density_Omega_G(arb_t res, const arb_t t_1, const arb_t s_1,
                                               void* param, const slong order, slong prec);

int interior_GW_current_energy_density_Omega_N_8(arb_t res, const arb_t t_1, const arb_t s_1,
                                                 const arb_t t_2, const arb_t s_2,
                                                 const arb_t t_3, const arb_t s_3,
                                                 const arb_t phi_12, const arb_t phi_23,
                                                 void* param, const slong order, slong prec);
int interior_GW_current_energy_density_Omega_P_8(arb_t res, const arb_t t_1, const arb_t s_1,
                                                 const arb_t t_2, const arb_t s_2,
                                                 const arb_t t_3, const arb_t s_3,
                                                 const arb_t phi_12, const arb_t phi_23,
                                                 void* param, const slong order, slong prec);
int interior_GW_current_energy_density_Omega_Z_5(arb_t res, const arb_t t_1, const arb_t s_1,
                                                 const arb_t t_2, const arb_t s_2, const arb_t phi_12,
                                                 void* param, const slong order, slong prec);
int interior_GW_current_energy_density_Omega_C_5(arb_t res, const arb_t t_1, const arb_t s_1,
                                                 const arb_t t_2, const arb_t s_2, const arb_t phi_12,
                                                 void* param, const slong order, slong prec);
int interior_GW_current_energy_density_Omega_R_6(arb_t res, const arb_t t_1, const arb_t s_1,
                                                 const arb_t t_2, const arb_t s_2,
                                                 const arb_t t_3, const arb_t s_3,
                                                 void* param, const slong order, slong prec);
int interior_GW_current_energy_density_Omega_H_4(arb_t res, const arb_t t_1, const arb_t s_1,
                                                 const arb_t t_2, const arb_t s_2,
                                                 void* param, const slong order, slong prec);


//方法三 2305.19950
int GW_current_energy_density_Omega_G(arb_t res, const arb_t k, slong prec); //当前的引力波能量密度

//下面是使用 cuba 中的积分方法来进行高维积分运算结果
int GW_current_energy_density_Omega_dim_8(arb_t O_N, arb_t O_P, const arb_t k, slong prec);//非高斯项
int GW_current_energy_density_Omega_dim_6(arb_t O_R, const arb_t k, slong prec);
int GW_current_energy_density_Omega_dim_5(arb_t O_Z, arb_t O_C, const arb_t k, slong prec);
int GW_current_energy_density_Omega_dim_4(arb_t O_H, const arb_t k, slong prec);

int GW_current_energy_density_Omega_dim_2(arb_t O_G, const arb_t k, slong prec);//高斯项

//当前的引力波能量密度，使用 cuba 积分的版本
//order = -1, 全部：高斯+非高斯
//order = 0, 高斯部分
//order = 2, 非高斯 (f_NL)^2 修正
//order = 4，非高斯 (f_NL)^4 修正
int GW_current_energy_density_cuba(arb_t res, const arb_t k, const int order, slong prec);


#endif // __PBHS_GW_INDUCED_THIRD_POWER_SPECTRA__  


