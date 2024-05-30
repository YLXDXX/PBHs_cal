#include "gw_induced.h"
#include <stdlib.h>

//注意，这里有两种方法，此两种方法给出的功率谱不相同
//但是当时间 η --> +∞ 时，此两种方法给出相同的功率谱 「Espinosa_01方法中需要设置辅助函数积分区间为[0,+∞]」
//故，对于当前的能量密度谱而言，是在 η --> +∞ 下给出的，此两种方法将给出相同的结果


//诱导引力波的功率谱
int GW_power_spectra(arb_t res, const arb_t eta, const arb_t k, slong prec)
{
    arb_t s;
    arb_init(s);
    
    switch(GW_induced_method)
    {
        case Espinosa_01 : //1804.07732
            GW_power_spectra_01(s, eta, k, prec);
            break;
        case Kohri_02 : //1804.08577
            GW_power_spectra_02(s, eta, k, prec);
            break;
        default :
            GW_power_spectra_01(s, eta, k, prec);
    }
    
    arb_set(res,s);
    
    arb_clear(s);
    return 0;
}

//当前的引力波能量密度
int GW_current_energy_density(arb_t res, const arb_t k, slong prec)
{
    arb_t s;
    arb_init(s);
    
    switch(GW_induced_method)
    {
        case Li_gauss :
            if( Int_GW_power_spectra_rectangle_adaptive )
            {
                //采用gauss_kronrod_iterate积分方法
                Integral_method_temp=Integral_method; //积分方法备份
                Integral_method=gauss_kronrod_iterate;
                
                GW_current_energy_density_Omega_G(s, k, prec);
                Integral_method=Integral_method_temp; //恢复之前的积分方法
                
            }else
            {
                GW_current_energy_density_Omega_G(s, k, prec);
            }
            
            break;
        case Espinosa_01 : //1804.07732
            GW_current_energy_density_01(s, k, prec);
            break;
        case Kohri_02 : //1804.08577
            GW_current_energy_density_02(s, k, prec);
            break;
        default :
            GW_current_energy_density_Omega_G(s, k, prec);
    }
    
    arb_set(res,s);
    
    arb_clear(s);
    return 0;
}


//当前的引力波能量密度，使用 cuba 积分的版本
int GW_current_energy_density_cuba(arb_t res, const arb_t k, const int order, slong prec)
{
    arb_t s,t,w,sum;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(sum);
    
    arb_zero(sum);//初始化为零
    
    switch( order )
    {
        case 0 : // order = 0, 高斯部分，只有一项
            
            GW_current_energy_density_Omega_dim_2(sum,k,prec);//Ω^G
            
            break;
        case 2 : // order = 2, 非高斯 (f_NL)^2 修正，三项
            
            GW_current_energy_density_Omega_dim_4(s,k,prec); //Ω^H
            GW_current_energy_density_Omega_dim_5(t,w,k,prec); //Ω^Z, Ω^C
            
            arb_add(sum,sum,s,prec);
            arb_add(sum,sum,t,prec);
            arb_add(sum,sum,w,prec);
            
            break;
        case 4 : // order = 4，非高斯 (f_NL)^4 修正，三项
            
            GW_current_energy_density_Omega_dim_6(s,k,prec); //Ω^R
            GW_current_energy_density_Omega_dim_8(t,w,k,prec); //Ω^N, Ω^P
            
            arb_add(sum,sum,s,prec);
            arb_add(sum,sum,t,prec);
            arb_add(sum,sum,w,prec);
            
            break;
        case -1 : //考虑所有：高斯+非高斯
            GW_current_energy_density_Omega_dim_2(s,k,prec);//Ω^G
            arb_add(sum,sum,s,prec);
            
            GW_current_energy_density_Omega_dim_4(s,k,prec); //Ω^H
            GW_current_energy_density_Omega_dim_5(t,w,k,prec); //Ω^Z, Ω^C
            arb_add(sum,sum,s,prec);
            arb_add(sum,sum,t,prec);
            arb_add(sum,sum,w,prec);
            
            GW_current_energy_density_Omega_dim_6(s,k,prec); //Ω^R
            GW_current_energy_density_Omega_dim_8(t,w,k,prec); //Ω^N, Ω^P
            arb_add(sum,sum,s,prec);
            arb_add(sum,sum,t,prec);
            arb_add(sum,sum,w,prec);
            
        default:
            printf("order 输入有误：\n\torder = -1, 全部：高斯+非高斯\n\torder = 0, 高斯部分\n\torder = 2, 非高斯 (f_NL)^2 修正\n\torder = 4，非高斯 (f_NL)^4 修正\n");
            exit(1);
    }
    
    arb_set(res,sum);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(sum);
    
    return 0;
}

