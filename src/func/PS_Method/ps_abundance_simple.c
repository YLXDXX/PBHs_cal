#include "ps_abundance_simple.h" 
#include <stdlib.h>

//使用Press-Schechter方法，进行简单的估算
//包括 密度扰动 δ，曲率扰动 ζ，compaction function C
//其中 δ 和 ζ 只能处理高斯的情况，而 C 可以处理非高斯情况
//在这里，关于不同视界质量的计算，统一都是利用窗口函数来获得

struct Window_func_and_R
{
    //利用R和窗口函数求方差时用
    arb_t R;
    enum WINDOW_FUNC_TYPE w_type;
    enum WINDOW_FUNC_TYPE w_type_second;
    int method;
};


//从PBHs的质量 M 得到相应的视界尺度 R 
void PS_variance_help_M_to_R(arb_t res, const arb_t m, slong prec)
{
    arb_t t,s,k;
    
    arb_init(t);
    arb_init(s);
    arb_init(k);
    
    //注意，为与compaction function 和peak theory保持一致
    //这里我们使用的质量，依然是相对质量 m=M/M_H 
    //其中的 M_H 取的是功率谱某一特征尺度 k_* 所对应的视界质量，对称谱，一般选为中心尺度 
    //参考δ谱的情况，k与R的转换关系为 R=2.743707269992269382559/k （单位Mpc）
    
    // k=K^{1/2}*[g(k_star)/g(k)]^{1/12}*k_star/m^{1/2} 
    //考虑到 [g(k_star)/g(k)]^{1/12} 大小在 1.2 左右，为计算的方面，我们忽略此因子
    //Func_k_to_degrees_of_freedom(t, s, k, prec);
    
    //简化为 k=K^{1/2}*k_star/m^{1/2}
    
    arb_sqrt(s,m,prec);
    arb_div(s,K_star,s,prec);
    
    arb_sqr(t,Mass_K,prec); //考虑了PBHs质量与视界质量的关系
    arb_mul(k,s,t,prec);
    
    arb_div(res, Delta_spectrum_x_m, k, prec);
    
    arb_clear(t);
    arb_clear(s);
    arb_clear(k);
}

//密度扰动的功率谱，以 ln(k) 作为变量
static void power_spectrum_density_contrast(arb_t res, const arb_t k, const arb_t R, slong prec)
{
    arb_t t,s,w;
    
    arb_init(t);
    arb_init(s);
    arb_init(w);
    
    //此处，对于密度扰动δ与曲率扰动ζ功率谱间的关系，使用简单的线性关系
    //对于w=1/3，P_δ(k)=16/81*(k*R)^4*P_ζ(k)
    
    arb_exp(w,k,prec); //w=e^k 恢复为原本的 k
    
    arb_mul(s,w,R,prec);
    arb_pow_ui(s,s,4,prec);
    arb_mul_ui(s,s,16,prec);
    arb_div_ui(s,s,81,prec);
    
    power_spectrum(t, k, prec); //以 ln(k) 作为变量
    
    arb_mul(res,s,t,prec);
    
    arb_clear(t);
    arb_clear(s);
    arb_clear(w);
}


static int interior_PS_variance_with_window_func(arb_t res, const arb_t k,
                                                 void* parameter, const slong order, slong prec)
{
    arb_t t,s,w;
    
    arb_init(t);
    arb_init(s);
    arb_init(w);
    
    struct Window_func_and_R *para; //这里不再需要分配内存
    para=parameter; //传入两个参数
    
    //注意，这里的积分变量k是取了对数的 ln(k)
    arb_exp(w,k,prec); //w=e^k 恢复为原本的 k，供窗口函数和转移函数调用
    
    if(order==0)
    {
        //方差 σ =∫W^2(k,R)*P(k)*T(k,η)*dln(k)
        
        Power_spectra_window_function_k(s, w, para->R, para->w_type, prec);
        arb_sqr(s,s,prec);
        
        if( para->method==0 ) //使用曲率扰动功率谱计算
        {
            power_spectrum(t, k, prec); //以 ln(k) 作变量
        }else //使用密度扰动的功率谱计算
        {
            power_spectrum_density_contrast(t, k, para->R, prec);
        }
        arb_mul(s,s,t,prec);
        
        //是否考虑转移函数的影响
        if( Transfer_Function )
        {
            Power_spectra_linear_transfer_function(t, w, para->R, prec); //传入的w已恢复成原本的k
            arb_sqr(t,t,prec);
            
            arb_mul(s,s,t,prec);
        }
        
        arb_set(res,s);
        
    }else //专门为compaction function Σ_{XY} 准备
    {
        //协方差 σ_{XY} =∫W_x(k,R)*W_y(k,R)*P(k)*T(k,η)*dln(k)
        
        Power_spectra_window_function_k(s, w, para->R, para->w_type, prec);
        Power_spectra_window_function_k(t, w, para->R, para->w_type_second, prec);
        arb_mul(s,s,t,prec);
        
        power_spectrum(t, k, prec); //以 ln(k) 作变量，为ζ的功率谱
        arb_mul(s,s,t,prec);
        
        //是否考虑转移函数的影响
        if( Transfer_Function )
        {
            Power_spectra_linear_transfer_function(t, w, para->R, prec); //传入的w已恢复成原本的k
            arb_sqr(t,t,prec);
            
            arb_mul(s,s,t,prec);
        }
        
        arb_set(res,s);
    }
    
    arb_clear(t);
    arb_clear(s);
    arb_clear(w);
    
    return 0;
}


//再利用窗口函数，得到质量 M 对应的统计量方差
void PS_variance_with_window_func(arb_t res, const arb_t R,
                                  const enum WINDOW_FUNC_TYPE w_type, const enum WINDOW_FUNC_TYPE w_type_second,
                                  const int method, slong prec)
{
    arb_t t,a,b;
    
    arb_init(t);
    arb_init(a);
    arb_init(b);
    
    int order;
    
    //这里，对于结构体 Window_func_and_R 需手动分配内存
    struct Window_func_and_R *para = (struct Window_func_and_R *)calloc(1,sizeof(struct Window_func_and_R));
    
    arb_init(para->R);//使用arb_t变量前初始化
    arb_set(para->R,R);
    
    para->method=method; // 0 表示默认使用曲率扰动ζ的功率谱计算，1 表示使用密度扰动的功率谱计算
    
    if(w_type==w_type_second)
    {
        para->w_type=w_type;
        order=0;
    }else
    {
        para->w_type=w_type;
        para->w_type_second=w_type_second;
        order=1;
    }
    
    //复用 compaction function 的协方差积分参数
    Integration_arb(t, interior_PS_variance_with_window_func, para, order, 
                    PS_Int_variance_min, PS_Int_variance_max, PS_Int_variance_precision,
                    Integration_iterate_min,Integration_iterate_max, prec);
    
    arb_set(res,t);
    
    arb_clear(t);
    arb_clear(a);
    arb_clear(b);
    free(para); //手动释放自定义结构体内存
}



//利用曲率扰动 ζ 估算，高斯情况


//利用密度扰动 δ 估算


//利用compaction function C 估算
