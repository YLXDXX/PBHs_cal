#include "basis.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cubaq.h>
#include <arb.h>
#include <quadmath.h> 

#define UNUSED(expr) do { (void)(expr); } while (0)

//利用 cuba 库进行多重积分的计算


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

/*
 * //有 ll 版本，即 long long int
 * void Vegas(const int ndim, const int ncomp,
 *           integrand_t integrand, void *userdata, const int nvec,
 *           const cubareal epsrel, const cubareal epsabs,
 *           const int flags, const int seed,
 *           const int mineval, const int maxeval,
 *           const int nstart, const int nincrease, const int nbatch,
 *           const int gridno, const char *statefile, void *spin,
 *           int *neval, int *fail,
 *           cubareal integral[], cubareal error[], cubareal prob[]);
 * 
 * void Suave(const int ndim, const int ncomp,
 *           integrand_t integrand, void *userdata, const int nvec,
 *           const cubareal epsrel, const cubareal epsabs,
 *           const int flags, const int seed,
 *           const int mineval, const int maxeval,
 *           const int nnew, const int nmin,
 *           const cubareal flatness, const char *statefile, void *spin,
 *           int *nregions, int *neval, int *fail,
 *           cubareal integral[], cubareal error[], cubareal prob[]);
 * 
 * void Divonne(const int ndim, const int ncomp,
 *             integrand_t integrand, void *userdata, const int nvec,
 *             const cubareal epsrel, const cubareal epsabs,
 *             const int flags, const int seed,
 *             const int mineval, const int maxeval,
 *             const int key1, const int key2, const int key3, const int maxpass,
 *             const cubareal border, const cubareal maxchisq, const cubareal mindeviation,
 *             const int ngiven, const int ldxgiven, cubareal xgiven[],
 *             const int nextra, peakfinder_t peakfinder,
 *             const char *statefile, void *spin,
 *             int *nregions, int *neval, int *fail,
 *             cubareal integral[], cubareal error[], cubareal prob[]);
 * 
 * void Cuhre(const int ndim, const int ncomp,
 *           integrand_t integrand, void *userdata, const int nvec,
 *           const cubareal epsrel, const cubareal epsabs,
 *           const int flags, const int mineval, const int maxeval,
 *           const int key,
 *           const char *statefile, void *spin,
 *           int *nregions, int *neval, int *fail,
 *           cubareal integral[], cubareal error[], cubareal prob[]);
 */

//
//输入参数
//
//ndim 维数
//ncomp 积分部件，即对于相同的区间，可以同时定义几个函数，一次计算几个积分
//integrand 被积函数，此函数定义为：
//    typedef int (*integrand_t)(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata);
//userdata 函数传递参数
//nvec 通常设为 1，与函数的矢量化相关
//epsrel, epsabs 相对误差和绝对误差，|I_c −I| ⩽ max(ε_abs, ε_rel|Ic|).
//flags 积分控制和日志打印
//seed 随机算法
//mineval, maxeval 最小最大计算次数
//statefile 中间态存储文件，可用于计算中断再恢复，一般设为""
//spin 子程序行为，一般设为 -1
//
//输出参数
//
//nregions 子区间数 Vegas算法没有该参数
//neval 积分计算实际次数
//fail 是否达到所需精度，达到反回 0
//integral 被积函数在单位超体积上的积分结果，[0,1]^dim
//error 对 integral 推测的绝对误差
//prob 对 integral的χ2-probability，可作为误差参考


//Vegas-specific Arguments
//nstart 每次迭代计算数
//nincrease 每次迭代被积函数积分量的增加
//nbatch 取样点批处理次数，一般设为 1000
//gridno 内部网格表量


static int Integrand_cuba(const int *ndim, const cubareal xx[],
                     const int *ncomp, cubareal ff[], void *userdata)
{
    //注意，Integration is always performed on the n-dimensional unit hypercube [O, 1]^n
    //需要变换区间，参见 https://giordano.github.io/Cuba.jl/stable/#Introduction
    
    //区间信息，通过 userdata 传递
    //U->interval_a;
    //U->interval_b;
    //U->func; //可有多个函数
    //*ncomp  //函数的个数
    //U->prec;
    
    struct User_Data *U;
    U=userdata;
    
    slong prec;
    prec=U->prec;
    
    arb_ptr coor,Jacobi; //各个维度坐标
    coor = _arb_vec_init(*ndim);
    Jacobi = _arb_vec_init(*ndim);
    
    arb_t t,s,a,b;
    
    arb_init(t);
    arb_init(s);
    arb_init(a);
    arb_init(b);
    
    int a_i,b_i;
    
    for( int i = 0; i < *ndim; i++ )
    {
        arb_set_d(coor+i, xx[i]); //将各个坐标转换到arb
        //对每个维度的区间做归一化
        arb_set(a,U->interval_a+i);
        arb_set(b,U->interval_b+i);
        
        a_i=arb_is_finite(a);
        b_i=arb_is_finite(b);
        
        
        //归一化区间
        if( a_i && b_i )
        {
            //[a,b]-->[0,1], x=a+(b−a)y, x --> a+(b−a)*x, Jacobi --> b-a
            arb_sub(t,b,a,prec);
            arb_set(Jacobi+i,t);
            arb_mul(t,t,coor+i,prec);
            arb_add(coor+i,t,a,prec);
        }else if ( a_i )
        {
            //[a,+∞] --> [0,1], x=a+y/(1−y), x --> a+x/(1-x), Jacobi --> 1/(1-x)^2
            arb_neg(t,coor+i);
            arb_add_ui(t,t,1,prec);
            arb_sqr(Jacobi+i,t,prec);
            arb_inv(Jacobi+i,Jacobi+i,prec);
            arb_div(t,coor+i,t,prec);
            arb_add(coor+i,a,t,prec);
        } else if ( b_i )
        {
            //[-∞,b] --> [0,1],  x=b+1-1/y, x --> b+1-1/x, Jacobi --> 1/x^2
            arb_inv(t,coor+i,prec);
            arb_sqr(Jacobi+i,t,prec);
            arb_neg(t,t);
            arb_add_ui(t,t,1,prec);
            arb_add(coor+i,b,t,prec);
        }else
        {
            //[-∞,+∞] --> [0,1], x --> (2*x-1)/[(1-x)*x], Jacobi (2*x^2-2*x+1)/[(1-x)^2 * x^2]
            arb_neg(t,coor+i);
            arb_add_ui(t,t,1,prec);
            arb_mul(t,t,coor+i,prec);
            
            arb_sqr(Jacobi,t,prec);
            arb_sqr(s,coor+i,prec);
            arb_sub(s,s,coor+i,prec);
            arb_mul_ui(s,s,2,prec);
            arb_add_ui(s,s,1,prec);
            arb_div(Jacobi,s,Jacobi,prec);
            
            arb_mul_ui(s,coor+i,2,prec);
            arb_sub_ui(s,s,1,prec);
            arb_add(coor+i,s,t,prec);
        }
    }
    
    for( int i = 0; i < *ncomp; i++ )
    {
        
        //维度不同，对应函数的参数个数不同
        switch (*ndim) 
        {
            case 8 :
                U->func[i](s,coor+0,coor+1,coor+2,coor+3,coor+4,coor+5,coor+6,coor+7,U->param,U->order,prec);
                break;
            case 7 :
                U->func[i](s,coor+0,coor+1,coor+2,coor+3,coor+4,coor+5,coor+6,U->param,U->order,prec);
                break;
            case 6 :
                U->func[i](s,coor+0,coor+1,coor+2,coor+3,coor+4,coor+5,U->param,U->order,prec);
                break;
            case 5 :
                U->func[i](s,coor+0,coor+1,coor+2,coor+3,coor+4,U->param,U->order,prec);
                break;
            case 4 :
                U->func[i](s,coor+0,coor+1,coor+2,coor+3,U->param,U->order,prec);
                break;
            case 3 :
                U->func[i](s,coor+0,coor+1,coor+2,U->param,U->order,prec);//三维函数
                break;
            case 2 :
                U->func[i](s,coor+0,coor+1,U->param,U->order,prec);
                break;
            case 1 :
                U->func[i](s,coor,U->param,U->order,prec);
                break;
            default:
                printf("未实现相应维数");
                exit(1);
        }
        
        //arb_printn(s, 50,0);printf("\n");
        
        for( int i = 0; i < *ndim; i++ ) //每个维度上，都要乘以相应的 Jacobi 因子
        {
            arb_mul(s,s,Jacobi+i,prec); 
        }
        
        ff[i]=strtoflt128( arb_get_str(s,35,ARB_STR_NO_RADIUS), 0 ); //将字符串转到 __float128
        //字符串转到 long double 使用 strtold()
        //printf("%i .... %.38Qg\n",i, ff[i] ); //打印 long double 用 %Lf, 打印 __float128 使用 %Qg
    }
    
    _arb_vec_clear(coor,*ndim);
    _arb_vec_clear(Jacobi,*ndim);
    
    arb_clear(t);
    arb_clear(s);
    arb_clear(a);
    arb_clear(b);
    
    return 0;
}



//设定积分区间
static void integration_interval_set(arb_ptr interval_a, arb_ptr interval_b,
                                     const char* t_upper,const int ndim, slong prec)
{
    char * two_pi="6.283185307";
    
    switch (ndim) 
    {
        case 8 :
            arb_set_str(interval_a+0,"0",prec);
            arb_set_str(interval_b+0,t_upper,prec);
            //arb_pos_inf(interval_b+0);
            
            arb_set_str(interval_a+1,"-1",prec);
            arb_set_str(interval_b+1,"1",prec);
            
            arb_set_str(interval_a+2,"0",prec);
            arb_set_str(interval_b+2,t_upper,prec);
            //arb_pos_inf(interval_b+2);
            
            arb_set_str(interval_a+3,"-1",prec);
            arb_set_str(interval_b+3,"1",prec);
            
            arb_set_str(interval_a+4,"0",prec);
            arb_set_str(interval_b+4,t_upper,prec);
            //arb_pos_inf(interval_b+4);
            
            arb_set_str(interval_a+5,"-1",prec);
            arb_set_str(interval_b+5,"1",prec);
            
            arb_set_str(interval_a+6,"0",prec);
            arb_set_str(interval_b+6,two_pi,prec); //2π
            
            arb_set_str(interval_a+7,"0",prec);
            arb_set_str(interval_b+7,two_pi,prec); //2π
            break;
        case 6 :
            arb_set_str(interval_a+0,"0",prec);
            arb_set_str(interval_b+0,t_upper,prec);
            
            arb_set_str(interval_a+1,"-1",prec);
            arb_set_str(interval_b+1,"1",prec);
            
            arb_set_str(interval_a+2,"0",prec);
            arb_set_str(interval_b+2,t_upper,prec);
            
            arb_set_str(interval_a+3,"-1",prec);
            arb_set_str(interval_b+3,"1",prec);
            
            arb_set_str(interval_a+4,"0",prec);
            arb_set_str(interval_b+4,t_upper,prec);
            
            arb_set_str(interval_a+5,"-1",prec);
            arb_set_str(interval_b+5,"1",prec);
            break;
        case 5 :
            arb_set_str(interval_a+0,"0",prec);
            arb_set_str(interval_b+0,t_upper,prec);
            
            arb_set_str(interval_a+1,"-1",prec);
            arb_set_str(interval_b+1,"1",prec);
            
            arb_set_str(interval_a+2,"0",prec);
            arb_set_str(interval_b+2,t_upper,prec);
            
            arb_set_str(interval_a+3,"-1",prec);
            arb_set_str(interval_b+3,"1",prec);
            
            arb_set_str(interval_a+4,"0",prec);
            arb_set_str(interval_b+4,two_pi,prec); //2π
            break;
        case 4 :
            arb_set_str(interval_a+0,"0",prec);
            arb_set_str(interval_b+0,t_upper,prec);
            
            arb_set_str(interval_a+1,"-1",prec);
            arb_set_str(interval_b+1,"1",prec);
            
            arb_set_str(interval_a+2,"0",prec);
            arb_set_str(interval_b+2,t_upper,prec);
            
            arb_set_str(interval_a+3,"-1",prec);
            arb_set_str(interval_b+3,"1",prec);
            break;
        case 3 :
            exit(0);
            break;
        case 2 :
            arb_set(interval_a+0,PS_Int_variance_min);
            arb_set(interval_b+0,PS_Int_variance_max);
            
            arb_set_str(interval_a+1,"0",prec);
            arb_set(interval_b+1,Pi);
            break;
        
        default:
            printf("维数错误，请检查\n");
            exit(1);
    }
}


static void res_print(int NCOMP,long long int neval, int fail,
                      cubareal integral[], cubareal error[], cubareal prob[])
{
    int comp;
    
    /*
    //对于积分结果，有相应的误差估计，加入球代数中
    char buf[256];
    char buf_zf[128];
    char *zf=" +/- ";
    quadmath_snprintf (buf, sizeof(buf), "%+-.15Qe", integral[0]);
    quadmath_snprintf (buf_zf, sizeof(buf_zf), "%.15Qe", error[0]);
    strcat(buf, zf);//连接字符串
    strcat(buf, buf_zf);
    */
    if( GW_dim_integral_res_print )
    {
        printf("VEGAS RESULT:\tneval %lld\tfail %d\n", neval, fail);
        for( comp = 0; comp < NCOMP; ++comp )
        {
            printf("VEGAS RESULT:\t%.12f +- %.12f\tp = %.12f\n",
                   (double)integral[comp], (double)error[comp], (double)prob[comp]);
        }
        printf("\n");
    }
}




static int interior_power_spectrum_non_Gaussian_f_Nl_3D(arb_t res, const arb_t p, const arb_t theta,
                                                 void* param, const slong order, slong prec)
{
    arb_t s,t,w,k,exp_p,exp_k,sqrt_kq;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(k);
    arb_init(exp_p);
    arb_init(exp_k);
    arb_init(sqrt_kq);
    
    //接收传入参数
    struct Func_transfer_parameter *func_k; //这里不需要手动分配，只需将其指向传入的指针即可
    func_k=param;
    
    arb_set(k,func_k->p_1); // k
    
    //注意到，这里传入的值，都是取对数后的值，包括 k 和 p
    //k=e^ln(k),  p=e^ln(p)
    
    arb_exp(exp_p,p,prec);
    arb_exp(exp_k,k,prec);
    
    
    //sqrt(k^2+p^2-2*k*p*cos(theta))
    arb_sqr(s,exp_k,prec);
    arb_sqr(t,exp_p,prec);
    arb_add(s,s,t,prec);
    
    arb_mul(w,exp_k,exp_p,prec);
    arb_mul_ui(w,w,2,prec);
    arb_cos(t,theta,prec);
    arb_mul(w,w,t,prec);
    
    arb_sub(s,s,w,prec);
    arb_sqrt(sqrt_kq,s,prec);
    
    
    //二元积分函数 开始
    
    //k^3/p^3*sin^2(θ)
    arb_pow_ui(s,exp_k,3,prec);
    arb_pow_ui(t,sqrt_kq,3,prec);
    arb_div(s,s,t,prec);
    arb_sin(t,theta,prec);
    arb_sqr(t,t,prec);
    arb_mul(s,s,t,prec);
    
    power_spectrum(t, p, prec); //高斯功率谱 P_ζ_G(k), 传入对数值
    arb_mul(s,s,t,prec);
    
    
    arb_log(t,sqrt_kq,prec);
    power_spectrum(w, t, prec); //高斯功率谱 P_ζ_G(k), 传入对数值
    
    
    arb_mul(res,s,w,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(k);
    arb_clear(exp_p);
    arb_clear(exp_k);
    arb_clear(sqrt_kq);
    
    return 0;
}


int power_spectrum_non_Gaussian_f_Nl(arb_t res, const arb_t k, slong prec)
{
    
    //cuba 积分参数控制
    /*********************************************************************/
    int NDIM = 2;
    int NCOMP = 1;
    void * USERDATA = NULL;
    int NVEC = 1;
    long double EPSREL = 1e-7; //相对精度
    long double EPSABS = 1e-12; //绝对精度
    int VERBOSE = 0;
    int LAST = 4;
    int SEED = 0;
    long int MINEVAL = 6000;
    long int MAXEVAL = 600000;
    
    long int NSTART = 1000;
    long int NINCREASE = 500;
    int NBATCH = 1000;
    int GRIDNO = 0;
    void * STATEFILE = NULL;
    void * SPIN = NULL;
    
    int NNEW = 1000;
    int NMIN = 2;
    long double FLATNESS = 25.;
    
    int KEY1 = 47;
    int KEY2 = 1;
    int KEY3 = 1;
    int MAXPASS = 5;
    long double BORDER = 0.;
    long double MAXCHISQ = 10.;
    long double MINDEVIATION = .25;
    int NGIVEN = 0;
    int LDXGIVEN = NDIM;
    int NEXTRA = 0;
    
    int KEY = 0;
    
    /*********************************************************************/
    UNUSED(KEY);UNUSED(NEXTRA);UNUSED(LDXGIVEN);UNUSED(NGIVEN);UNUSED(MINDEVIATION);UNUSED(MAXCHISQ);
    UNUSED(BORDER);UNUSED(MAXPASS);UNUSED(KEY3);UNUSED(KEY2);UNUSED(KEY1);UNUSED(FLATNESS);
    UNUSED(NMIN);UNUSED(NNEW);UNUSED(SPIN);UNUSED(STATEFILE);UNUSED(GRIDNO);UNUSED(NBATCH);
    UNUSED(NSTART);UNUSED(NINCREASE);UNUSED(LAST);UNUSED(USERDATA);UNUSED(SEED);
    
    arb_t s,t,w;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    struct User_Data *U;
    //分配内存
    U = (struct User_Data *)calloc(1,sizeof(struct User_Data));
    //初始化其中变量
    U->prec=prec;
    U->interval_a=_arb_vec_init(NDIM);
    U->interval_b=_arb_vec_init(NDIM);
    U->func = ( func_NULL * )calloc(NCOMP,sizeof( func_NULL )); //根据函数个数，初始化相应的函数指针数组
    U->func[0]=interior_power_spectrum_non_Gaussian_f_Nl_3D; //两个被积函数
    //U->func[1]=interior_GW_current_energy_density_Omega_P_8;
    
    //设定积分区间
    integration_interval_set(U->interval_a, U->interval_b, "1E3", NDIM, prec);
    
    //其中η和k是积分函数的参数，需传入
    //这里，对于结构体 Func_transfer_parameter 需手动分配内存
    struct Func_transfer_parameter *func_k;
    func_k=(struct Func_transfer_parameter *)calloc(1,sizeof(struct Func_transfer_parameter));
    arb_init(func_k->p_1);//使用arb_t变量前初始化
    
    arb_set(func_k->p_1,k); //设定传递参数
    U->param=func_k; //相应参数通过 User_Data 交互传入
    U->order=0;
    
    
    
    //积分前面的系数 72/25*π^3*(f_nl)^2
    arb_set_ui(s,72);
    arb_div_ui(s,s,25,prec);
    arb_pow_ui(t,Pi,3,prec);
    arb_mul(s,s,t,prec);
    arb_sqr(t,Power_expansion_f,prec); 
    arb_mul(s,s,t,prec);
    
    
    int nregions, fail;
    long long int neval;
    cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];
    
    //对于有限区间，采用 llDivonne, llCuhre 效果较好
    //对于无限区间，采用 llVegas, llSuave, llCuhre 效果较好
    
    /*
    llSuave(NDIM, NCOMP, Integrand_cuba, U, NVEC,
            EPSREL, EPSABS, VERBOSE | LAST, SEED,
            MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
            STATEFILE, SPIN,
            &nregions, &neval, &fail, integral, error, prob);
    */
    /*
    llDivonne(NDIM, NCOMP, Integrand_cuba, U, NVEC,
              EPSREL, EPSABS, VERBOSE, SEED,
              MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
              BORDER, MAXCHISQ, MINDEVIATION,
              NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
              STATEFILE, SPIN,
              &nregions, &neval, &fail, integral, error, prob);
    */
    
    llCuhre(NDIM, NCOMP, Integrand_cuba, U, NVEC,
            EPSREL, EPSABS, VERBOSE | LAST,
            MINEVAL, MAXEVAL, KEY,
            STATEFILE, SPIN,
            &nregions, &neval, &fail, integral, error, prob);
    
    res_print(NCOMP, neval, fail, integral, error, prob);//打印结果
    
    arb_set_d(t,integral[0]);//转入arb
    //arb_set_d(w,integral[1]);//转入arb
    //arb_set_str(t,buf,prec);
    
    arb_mul(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    
    arb_clear(func_k->p_1);
    free(func_k);
    
    return 0;
}



