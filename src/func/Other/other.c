#include "other.h" 
#include <stdlib.h>
#include <string.h> 

//Heaviside Theta function
int Heaviside_Theta_function(arb_t res,const arb_t x,slong prec)
{
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    
    arb_zero(t);
    
    if ( arb_gt(x,t) ) //x>0, f(x)=1
    {
        arb_one(res);
    }else if ( arb_lt(x,t) ) //x<0, f(x)=0
    {
        arb_zero(res);
    }else if ( arb_eq(x,t) )  //x=0, f(x)=1/2
    {
        arb_one(s); //1/2
        arb_div_ui(s,s,2,prec);
        
        arb_set(res,s);
    }else
    {
        printf("x 取值有误\n");
        exit(1);
    }
    
    
    arb_clear(s);
    arb_clear(t);
    
    return 0;
}


//诱导引力波的频率 f nHz 转成 波数 k Mpc^-1
void Func_GW_f_to_k(arb_t k, const arb_t f, slong prec)
{
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    //f=1.54*(k/10^6)
    arb_mul_ui(s,f,1E6, prec);
    arb_set_str(t,"1.54",prec);
    arb_div(k,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
}

