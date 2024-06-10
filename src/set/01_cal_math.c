#include "header/01_cal_math.h" 

void Set_cal_math(slong prec)
{
    //设定常用数学常数，供后面调用
    arb_const_pi(Pi, prec);//常数π
    arb_mul_si(Pi_2,Pi,2,prec);//常数2π (πx2)
    
    
    //找根行为设定
    Find_root_method=Brent_method; // Brent_method/bisection_methode
    
    //积分行为设定
    Integral_method=double_exponential; // gauss_kronrod_iterate/double_exponential
    
    Multithreaded_divide_integration_interval_number=16; //多线程计算区间分隔数目
    Multithreaded_number=11; //多线程计算线程总数
    
    
    arb_init(INT_MIN_INTERVAL); //积分最小子区间间隔，一般十进制设为 10^(-prec/4) ，主要为防止有些积分区间无限细分下去
    arb_set_str(INT_MIN_INTERVAL,"10",prec);
    long int t_int_gap=prec/4.5; //涉汲到整数除法，这里会自动取整，这里是向零取整，对正数而言为向下取整
    arb_pow_ui(INT_MIN_INTERVAL,INT_MIN_INTERVAL,t_int_gap,prec);
    arb_inv(INT_MIN_INTERVAL,INT_MIN_INTERVAL,prec);
    
    
    switch(Integral_method) //对不同的积分行为设定不同的积分参数
    {
        case gauss_kronrod_iterate : //gauss_kronrod迭代版本
            
            //2^12=4096 2^13=8192 2^15=32768 2^17=131072
            Integration_iterate_min=32; //积分最小迭代区间 interval_min = (b-a)/step_min
            Integration_iterate_max=25000; //积分最大迭代次数
            
            break;
            
        case double_exponential :
            Integration_iterate_min=4; //最少迭代次数
            Integration_iterate_max=16; //最大迭代次数
            break;
            
        default :
            exit(1);
    }
} 
