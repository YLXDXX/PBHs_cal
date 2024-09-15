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



//线性分割，类似于 numpy.linspace
void Get_interval_linspace_point(arb_ptr x, const arb_t a, const arb_t b, const slong N, slong prec)
{
    //在区间[a,b]中产生 N 个点（包含端点）
    //N个点在 x 上均匀分酆
    //[a,b]上取 N 个点，共分成 N-1 份
    
    arb_t s,delta;
    arb_init(s);
    arb_init(delta);
    
    arb_sub(delta,b,a,prec);
    arb_div_si(delta,delta,N-1,prec);
    
    for(slong i=0; i < N; i++ )
    {
        arb_mul_si(s,delta,i,prec);
        arb_add(x+i,a,s,prec);
    }
    
    arb_clear(s);
    arb_clear(delta);
}

//对数分割，类似于 numpy.logspace
void Get_interval_logspace_point(arb_ptr x, const arb_t a, const arb_t b, const slong N, slong prec)
{
    //很多图是对数图，其函数取点应是对数图log(x) 上均匀，但在坐标值 x 上不均匀
    //在区间[a,b]中产生 N 个点（包含端点）
    //[a,b]上取 N 个点，共分成 N-1 份
    //矢量 arb_ptr的维数致少为 N
    
    if( (!arb_is_positive(a)) || (!arb_is_positive(b)) )
    {
        printf("区间端点的符号有问题：\na=");
        arb_printn(a, 50,0);printf("\n");
        exit(1);
    }
        
    arb_t s,c_i,log_a,log_b,ln_10,log_delta;
    arb_init(s);
    arb_init(c_i);
    arb_init(log_a);
    arb_init(log_b);
    arb_init(ln_10);
    arb_init(log_delta);
    
    arb_set_si(s,10);
    arb_log(ln_10,s,prec);
    
    arb_log(s,a,prec);
    arb_div(log_a,s,ln_10,prec);
    
    arb_log(s,b,prec);
    arb_div(log_b,s,ln_10,prec);
    
    arb_sub(log_delta,log_b,log_a,prec); //对数上，各点间隔
    arb_div_si(log_delta,log_delta,N-1,prec);
    
    for(slong i=0; i < N; i++ )
    {
        //c_i=log(a) + [ log(b)-log(a) ]/(N-1)*i
        arb_mul_si(s,log_delta,i,prec);
        arb_add(c_i,log_a,s,prec);
        
        //x_i=Exp(c_i*ln(10))
        arb_mul(s,c_i,ln_10,prec);
        arb_exp(x+i,s,prec);
    }
    
    arb_clear(s);
    arb_clear(c_i);
    arb_clear(log_a);
    arb_clear(log_b);
    arb_clear(ln_10);
    arb_clear(log_delta);
}
