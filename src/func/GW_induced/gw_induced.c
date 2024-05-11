#include "gw_induced.h"

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
            GW_power_spectra_02(s, eta, k, prec);
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
        case Espinosa_01 : //1804.07732
            GW_current_energy_density_01(s, k, prec);
            break;
        case Kohri_02 : //1804.08577
            GW_current_energy_density_02(s, k, prec);
            break;
        default :
            GW_current_energy_density_02(s, k, prec);
    }
    
    arb_set(res,s);
    
    arb_clear(s);
    return 0;
}


