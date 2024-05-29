#include "cal_func.h" 
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


arb_ptr Fitting_a_i,Fitting_b_i,Fitting_c_i,Fitting_d_i; //自由度数拟合系数
arb_t Fitting_m_e,Fitting_m_mu,Fitting_m_pi_0,Fitting_m_pi_pm,
      Fitting_m_1,Fitting_m_2,Fitting_m_3,Fitting_m_4;//自由度数拟合系数
arb_ptr Fitting_func_f_rho,Fitting_func_b_rho,Fitting_func_f_s,Fitting_func_b_s,Fitting_func_S_fit;//自由度数拟合函数系数
arb_ptr Fitting_func_g_star_r,Fitting_func_g_star_s;
static int get_effective_degrees_of_freedom_fitting_coefficient(slong prec)
{
    if(Fitting_a_i != NULL)
    {   
        return 0;
    }
    //见 1803.01038 Table 1
    //自由度数拟合系数
    
    Fitting_a_i = _arb_vec_init(12);
    Fitting_b_i = _arb_vec_init(12);
    Fitting_c_i = _arb_vec_init(12);
    Fitting_d_i = _arb_vec_init(12);
    
    //拟全系数赋值
    arb_set_str(Fitting_a_i+0,"1",prec);
    arb_set_str(Fitting_a_i+1,"1.11724E+00",prec);
    arb_set_str(Fitting_a_i+2,"3.12672E-01",prec);
    arb_set_str(Fitting_a_i+3,"-4.68049E-02",prec);
    arb_set_str(Fitting_a_i+4,"-2.65004E-02",prec);
    arb_set_str(Fitting_a_i+5,"-1.19760E-03",prec);
    arb_set_str(Fitting_a_i+6,"1.82812E-04",prec);
    arb_set_str(Fitting_a_i+7,"1.36436E-04",prec);
    arb_set_str(Fitting_a_i+8,"8.55051E-05",prec);
    arb_set_str(Fitting_a_i+9,"1.22840E-05",prec);
    arb_set_str(Fitting_a_i+10,"3.82259E-07",prec);
    arb_set_str(Fitting_a_i+11,"-6.87035E-09",prec);
    
    arb_set_str(Fitting_b_i+0,"1.43382E-02",prec);
    arb_set_str(Fitting_b_i+1,"1.37559E-02",prec);
    arb_set_str(Fitting_b_i+2,"2.92108E-03",prec);
    arb_set_str(Fitting_b_i+3,"-5.38533E-04",prec);
    arb_set_str(Fitting_b_i+4,"-1.62496E-04",prec);
    arb_set_str(Fitting_b_i+5,"-2.87906E-05",prec);
    arb_set_str(Fitting_b_i+6,"-3.84278E-06",prec);
    arb_set_str(Fitting_b_i+7,"2.78776E-06",prec);
    arb_set_str(Fitting_b_i+8,"7.40342E-07",prec);
    arb_set_str(Fitting_b_i+9,"1.17210E-07",prec);
    arb_set_str(Fitting_b_i+10,"3.72499E-09",prec);
    arb_set_str(Fitting_b_i+11,"-6.74107E-11",prec);
    
    arb_set_str(Fitting_c_i+0,"1",prec);
    arb_set_str(Fitting_c_i+1,"6.07869E-01",prec);
    arb_set_str(Fitting_c_i+2,"-1.54485E-01",prec);
    arb_set_str(Fitting_c_i+3,"-2.24034E-01",prec);
    arb_set_str(Fitting_c_i+4,"-2.82147E-02",prec);
    arb_set_str(Fitting_c_i+5,"2.90620E-02",prec);
    arb_set_str(Fitting_c_i+6,"6.86778E-03",prec);
    arb_set_str(Fitting_c_i+7,"-1.00005E-03",prec);
    arb_set_str(Fitting_c_i+8,"-1.69104E-04",prec);
    arb_set_str(Fitting_c_i+9,"1.06301E-05",prec);
    arb_set_str(Fitting_c_i+10,"1.69528E-06",prec);
    arb_set_str(Fitting_c_i+11,"-9.33311E-08",prec);
    
    arb_set_str(Fitting_d_i+0,"7.07388E+01",prec);
    arb_set_str(Fitting_d_i+1,"9.18011E+01",prec);
    arb_set_str(Fitting_d_i+2,"3.31892E+01",prec);
    arb_set_str(Fitting_d_i+3,"-1.39779E+00",prec);
    arb_set_str(Fitting_d_i+4,"-1.52558E+00",prec);
    arb_set_str(Fitting_d_i+5,"-1.97857E-02",prec);
    arb_set_str(Fitting_d_i+6,"-1.60146E-01",prec);
    arb_set_str(Fitting_d_i+7,"8.22615E-05",prec);
    arb_set_str(Fitting_d_i+8,"2.02651E-02",prec);
    arb_set_str(Fitting_d_i+9,"-1.82134E-05",prec);
    arb_set_str(Fitting_d_i+10,"7.83943E-05",prec);
    arb_set_str(Fitting_d_i+11,"7.13518E-05",prec);
    
    //对应质量
    arb_set_str(Fitting_m_e,"511E-6",prec);
    arb_set_str(Fitting_m_mu,"0.1056",prec);
    arb_set_str(Fitting_m_pi_0,"0.135",prec);
    arb_set_str(Fitting_m_pi_pm,"0.140",prec);
    arb_set_str(Fitting_m_1,"0.5",prec);
    arb_set_str(Fitting_m_2,"0.77",prec);
    arb_set_str(Fitting_m_3,"1.2",prec);
    arb_set_str(Fitting_m_4,"2",prec);
    
    //拟合函数系数
    Fitting_func_f_rho = _arb_vec_init(4);
    Fitting_func_b_rho = _arb_vec_init(4);
    Fitting_func_f_s = _arb_vec_init(4);
    Fitting_func_b_s = _arb_vec_init(4);
    Fitting_func_S_fit = _arb_vec_init(4);
    
    //函数的系数从左到右依次4个实数
    arb_set_str(Fitting_func_f_rho+0,"-1.04855",prec);//f_rho(x)
    arb_set_str(Fitting_func_f_rho+1,"1.03757",prec);
    arb_set_str(Fitting_func_f_rho+2,"0.508630",prec);
    arb_set_str(Fitting_func_f_rho+3,"0.0893988",prec);
    
    arb_set_str(Fitting_func_b_rho+0,"-1.03149",prec);//b_rho(x)
    arb_set_str(Fitting_func_b_rho+1,"1.03317",prec);
    arb_set_str(Fitting_func_b_rho+2,"0.398264",prec);
    arb_set_str(Fitting_func_b_rho+3,"0.0648056",prec);
    
    arb_set_str(Fitting_func_f_s+0,"-1.04190",prec);//f_s(x)
    arb_set_str(Fitting_func_f_s+1,"1.03400",prec);
    arb_set_str(Fitting_func_f_s+2,"0.456426",prec);
    arb_set_str(Fitting_func_f_s+3,"0.0595248",prec);
    
    arb_set_str(Fitting_func_b_s+0,"-1.03365",prec);//b_s(x)
    arb_set_str(Fitting_func_b_s+1,"1.03397",prec);
    arb_set_str(Fitting_func_b_s+2,"0.342548",prec);
    arb_set_str(Fitting_func_b_s+3,"0.0506182",prec);
    
    arb_set_str(Fitting_func_S_fit+0,"-1.0419",prec);//S_fit(x)
    arb_set_str(Fitting_func_S_fit+1,"1.034",prec);
    arb_set_str(Fitting_func_S_fit+2,"0.456426",prec);
    arb_set_str(Fitting_func_S_fit+3,"0.0595249",prec);
    
    //自由度数函数系数
    Fitting_func_g_star_r=_arb_vec_init(10);
    Fitting_func_g_star_s=_arb_vec_init(10);
    
    arb_set_str(Fitting_func_g_star_r+0,"2.030",prec);//g_{*,ρ}(x)
    arb_set_str(Fitting_func_g_star_r+1,"1.353",prec);
    arb_set_str(Fitting_func_g_star_r+2,"3.495",prec);
    arb_set_str(Fitting_func_g_star_r+3,"3.446",prec);
    arb_set_str(Fitting_func_g_star_r+4,"1.05",prec);
    arb_set_str(Fitting_func_g_star_r+5,"2.08",prec);
    arb_set_str(Fitting_func_g_star_r+6,"4.165",prec);
    arb_set_str(Fitting_func_g_star_r+7,"30.55",prec);
    arb_set_str(Fitting_func_g_star_r+8,"89.4",prec);
    arb_set_str(Fitting_func_g_star_r+9,"8209",prec);
    
    arb_set_str(Fitting_func_g_star_s+0,"2.008",prec);//g_{*,s}(x)
    arb_set_str(Fitting_func_g_star_s+1,"1.923",prec);
    arb_set_str(Fitting_func_g_star_s+2,"3.442",prec);
    arb_set_str(Fitting_func_g_star_s+3,"3.468",prec);
    arb_set_str(Fitting_func_g_star_s+4,"1.034",prec);
    arb_set_str(Fitting_func_g_star_s+5,"2.068",prec);
    arb_set_str(Fitting_func_g_star_s+6,"4.16",prec);
    arb_set_str(Fitting_func_g_star_s+7,"30.55",prec);
    arb_set_str(Fitting_func_g_star_s+8,"90",prec);
    arb_set_str(Fitting_func_g_star_s+9,"6209",prec);
    
    return 0;
}

//拟合所用函数
static void f_rho(arb_t res, const arb_t x, slong prec)
{
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    //右边括号
    arb_mul(t,Fitting_func_f_rho+1,x,prec);
    arb_add_ui(t,t,1,prec);
    
    arb_sqr(s,x,prec);
    arb_mul(s,Fitting_func_f_rho+2,s,prec);
    arb_add(t,t,s,prec);
    
    arb_pow_ui(s,x,3,prec);
    arb_mul(s,Fitting_func_f_rho+3,s,prec);
    arb_add(t,t,s,prec);
    
    //左边指数
    arb_mul(s,Fitting_func_f_rho+0,x,prec);
    arb_exp(s,s,prec);
    
    arb_mul(res,t,s,prec);
    
    arb_clear(s);
    arb_clear(t);
}

static void b_rho(arb_t res, const arb_t x, slong prec)
{
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    //右边括号
    arb_mul(t,Fitting_func_b_rho+1,x,prec);
    arb_add_ui(t,t,1,prec);
    
    arb_sqr(s,x,prec);
    arb_mul(s,Fitting_func_b_rho+2,s,prec);
    arb_add(t,t,s,prec);
    
    arb_pow_ui(s,x,3,prec);
    arb_mul(s,Fitting_func_b_rho+3,s,prec);
    arb_add(t,t,s,prec);
    
    //左边指数
    arb_mul(s,Fitting_func_b_rho+0,x,prec);
    arb_exp(s,s,prec);
    
    arb_mul(res,t,s,prec);
    
    arb_clear(s);
    arb_clear(t);
}

static void f_s(arb_t res, const arb_t x, slong prec)
{
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    //右边括号
    arb_mul(t,Fitting_func_f_s+1,x,prec);
    arb_add_ui(t,t,1,prec);
    
    arb_sqr(s,x,prec);
    arb_mul(s,Fitting_func_f_s+2,s,prec);
    arb_add(t,t,s,prec);
    
    arb_pow_ui(s,x,3,prec);
    arb_mul(s,Fitting_func_f_s+3,s,prec);
    arb_add(t,t,s,prec);
    
    //左边指数
    arb_mul(s,Fitting_func_f_s+0,x,prec);
    arb_exp(s,s,prec);
    
    arb_mul(res,t,s,prec);
    
    arb_clear(s);
    arb_clear(t);
}

static void b_s(arb_t res, const arb_t x, slong prec)
{
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    //右边括号
    arb_mul(t,Fitting_func_b_s+1,x,prec);
    arb_add_ui(t,t,1,prec);
    
    arb_sqr(s,x,prec);
    arb_mul(s,Fitting_func_b_s+2,s,prec);
    arb_add(t,t,s,prec);
    
    arb_pow_ui(s,x,3,prec);
    arb_mul(s,Fitting_func_b_s+3,s,prec);
    arb_add(t,t,s,prec);
    
    //左边指数
    arb_mul(s,Fitting_func_b_s+0,x,prec);
    arb_exp(s,s,prec);
    
    arb_mul(res,t,s,prec);
    
    arb_clear(s);
    arb_clear(t);
}

static void S_fit(arb_t res, const arb_t x, slong prec)
{
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    //右边括号
    arb_mul(t,Fitting_func_S_fit+1,x,prec);
    arb_add_ui(t,t,1,prec);
    
    arb_sqr(s,x,prec);
    arb_mul(s,Fitting_func_S_fit+2,s,prec);
    arb_add(t,t,s,prec);
    
    arb_pow_ui(s,x,3,prec);
    arb_mul(s,Fitting_func_S_fit+3,s,prec);
    arb_add(t,t,s,prec);
    
    //左边指数
    arb_mul(s,Fitting_func_S_fit+0,x,prec);
    arb_exp(s,s,prec);
    
    arb_mul(t,t,s,prec);
    
    arb_mul_ui(t,t,7,prec);
    arb_div_ui(t,t,4,prec);
    
    arb_add_ui(res,t,1,prec);
    
    arb_clear(s);
    arb_clear(t);
}

//计算对应温度的相对论自由度数和熵自由度数
//Fitting functions for effective degrees of freedom
//详情可参考 1803.01038
void Effective_degrees_of_freedom_fit(arb_t res_r, arb_t res_s, const arb_t T, char * unit, slong prec)
{
    get_effective_degrees_of_freedom_fitting_coefficient(prec); //获取拟合系数
    
    arb_t t,tt,TT,s,w,sum_1,sum_2,sum_3,sum_4;
    arb_init(t);
    arb_init(tt);
    arb_init(TT);
    arb_init(s);
    arb_init(w);
    arb_init(sum_1);
    arb_init(sum_2);
    arb_init(sum_3);
    arb_init(sum_4);
    
    arb_zero(sum_1);
    arb_zero(sum_2);
    arb_zero(sum_3);
    arb_zero(sum_4);
    
    
    //单位统一为 Gev
    if( !strcmp(unit,"Gev") )
    {
        arb_set(TT,T);
    }else if ( !strcmp(unit,"Mev") )
    {
        arb_div_ui(TT,T,1000,prec); //Mev --> Gev
    }else
    {
        printf("单位输入错误，应为 Mev 或 Gev\n");
        exit(1);
    }
    
    
    long double temp;
    temp=0.120; //0.120Gev
    arb_set_d(t,temp);
    if( arb_gt(TT,t) ) //拟合范围 120 MeV ≤T ≤10^16 GeV
    {
        //t= ln(T[GeV]).
        arb_log(t,TT,prec); //这里是以 Gev 10^9 为单位；Mev 10^6
        
        for( int i = 0; i < 12; i++ )
        {
            arb_pow_ui(tt,t,i,prec);
            arb_mul(s,Fitting_a_i+i,tt,prec);
            arb_mul(w,Fitting_b_i+i,tt,prec);
            
            arb_add(sum_1,sum_1,s,prec);//分子
            arb_add(sum_2,sum_2,w,prec);//分母
            
            arb_mul(s,Fitting_c_i+i,tt,prec);
            arb_mul(w,Fitting_d_i+i,tt,prec);
            
            arb_add(sum_3,sum_3,s,prec);//分子
            arb_add(sum_4,sum_4,w,prec);//分母
        }
        
        arb_div(res_r,sum_1,sum_2,prec);//相对论自由度数
        
        arb_div(s,sum_3,sum_4,prec);
        arb_add_ui(s,s,1,prec);
        arb_div(res_s,res_r,s,prec); //熵自由度数
    }else
    {
        //相对论自由度数
        arb_div(t,Fitting_m_e,TT,prec);
        S_fit(s,t,prec);
        
        arb_one(w);
        arb_mul_ui(w,w,4,prec);
        arb_div_ui(w,w,3,prec);
        arb_pow(s,s,w,prec);
        
        arb_mul(t,Fitting_func_g_star_r+1,s,prec);
        arb_add(sum_1,t,Fitting_func_g_star_r+0,prec);
        
        arb_div(t,Fitting_m_e,TT,prec);
        f_rho(s,t,prec);
        arb_mul(t,Fitting_func_g_star_r+2,s,prec);
        arb_add(sum_1,sum_1,t,prec);
        
        arb_div(t,Fitting_m_mu,TT,prec);
        f_rho(s,t,prec);
        arb_mul(t,Fitting_func_g_star_r+3,s,prec);
        arb_add(sum_1,sum_1,t,prec);
        
        arb_div(t,Fitting_m_pi_0,TT,prec);
        b_rho(s,t,prec);
        arb_mul(t,Fitting_func_g_star_r+4,s,prec);
        arb_add(sum_1,sum_1,t,prec);
        
        arb_div(t,Fitting_m_pi_pm,TT,prec);
        b_rho(s,t,prec);
        arb_mul(t,Fitting_func_g_star_r+5,s,prec);
        arb_add(sum_1,sum_1,t,prec);
        
        arb_div(t,Fitting_m_1,TT,prec);
        b_rho(s,t,prec);
        arb_mul(t,Fitting_func_g_star_r+6,s,prec);
        arb_add(sum_1,sum_1,t,prec);
        
        arb_div(t,Fitting_m_2,TT,prec);
        b_rho(s,t,prec);
        arb_mul(t,Fitting_func_g_star_r+7,s,prec);
        arb_add(sum_1,sum_1,t,prec);
        
        arb_div(t,Fitting_m_3,TT,prec);
        b_rho(s,t,prec);
        arb_mul(t,Fitting_func_g_star_r+8,s,prec);
        arb_add(sum_1,sum_1,t,prec);
        
        arb_div(t,Fitting_m_4,TT,prec);
        b_rho(s,t,prec);
        arb_mul(t,Fitting_func_g_star_r+9,s,prec);
        arb_add(sum_1,sum_1,t,prec);
        
        arb_set(res_r,sum_1);
        
        //熵自由度数
        arb_div(t,Fitting_m_e,TT,prec);
        S_fit(s,t,prec);
        arb_mul(t,Fitting_func_g_star_s+1,s,prec);
        arb_add(sum_2,t,Fitting_func_g_star_s+0,prec);
        
        arb_div(t,Fitting_m_e,TT,prec);
        f_s(s,t,prec);
        arb_mul(t,Fitting_func_g_star_s+2,s,prec);
        arb_add(sum_2,sum_2,t,prec);
        
        arb_div(t,Fitting_m_mu,TT,prec);
        f_s(s,t,prec);
        arb_mul(t,Fitting_func_g_star_s+3,s,prec);
        arb_add(sum_2,sum_2,t,prec);
        
        arb_div(t,Fitting_m_pi_0,TT,prec);
        b_s(s,t,prec);
        arb_mul(t,Fitting_func_g_star_s+4,s,prec);
        arb_add(sum_2,sum_2,t,prec);
        
        arb_div(t,Fitting_m_pi_pm,TT,prec);
        b_s(s,t,prec);
        arb_mul(t,Fitting_func_g_star_s+5,s,prec);
        arb_add(sum_2,sum_2,t,prec);
        
        arb_div(t,Fitting_m_1,TT,prec);
        b_s(s,t,prec);
        arb_mul(t,Fitting_func_g_star_s+6,s,prec);
        arb_add(sum_2,sum_2,t,prec);
        
        arb_div(t,Fitting_m_2,TT,prec);
        b_s(s,t,prec);
        arb_mul(t,Fitting_func_g_star_s+7,s,prec);
        arb_add(sum_2,sum_2,t,prec);
        
        arb_div(t,Fitting_m_3,TT,prec);
        b_s(s,t,prec);
        arb_mul(t,Fitting_func_g_star_s+8,s,prec);
        arb_add(sum_2,sum_2,t,prec);
        
        arb_div(t,Fitting_m_4,TT,prec);
        b_s(s,t,prec);
        arb_mul(t,Fitting_func_g_star_s+9,s,prec);
        arb_add(sum_2,sum_2,t,prec);
        
        arb_set(res_s,sum_2);
    }
    
    
    arb_clear(t);
    arb_clear(tt);
    arb_clear(TT);
    arb_clear(s);
    arb_clear(w);
    arb_clear(sum_1);
    arb_clear(sum_2);
    arb_clear(sum_3);
    arb_clear(sum_4);
}


//求根函数
static int interior_k_to_degrees_of_freedom(arb_t res, const arb_t T, void *k, const slong order, slong prec)
{
    //利用 k 与 温度 T 的关系, 详细公式见 1812.00674 (7)
    //通过 k 找到对应的温度，再通过温度 T 求出对应的自由度数
    //为了减小温度的求根区间，这里传入的温度是 ln(T), [-12, +12]
    
    arb_t T_ratio,k_ratio,g,g_s,t,s,w;
    arb_init(T_ratio);
    arb_init(k_ratio);
    arb_init(g);
    arb_init(g_s);
    arb_init(t);
    arb_init(s);
    arb_init(w);
    
    
    arb_exp(w,T,prec); //温度转换
    arb_div(T_ratio,w,T_scale_eq,prec);
    
    arb_div(k_ratio,k,K_scale_eq,prec); //这里单位为 Mpc^-1
    
    arb_div(s,T_ratio,k_ratio,prec); //通过比值，避免求根中大数的产生
    
    Effective_degrees_of_freedom_fit(g, g_s, w,"Gev",prec);
    
    arb_one(t);
    arb_div_ui(t,t,3,prec);
    arb_div(w,effective_g_star_eq_entropy,g_s,prec);
    arb_pow(w,w,t,prec);
    arb_mul(s,s,w,prec);
    
    arb_one(t);
    arb_div_ui(t,t,2,prec);
    arb_div(w,g,effective_g_star_eq,prec);
    arb_pow(w,w,t,prec);
    arb_mul(s,s,w,prec);
    
    arb_one(t); //最前面的系数
    arb_mul_ui(t,t,2,prec);
    arb_sqrt(t,t,prec);
    arb_sub_ui(t,t,1,prec);
    arb_mul_ui(t,t,2,prec);
    arb_mul(s,s,t,prec);
    
    arb_sub_ui(res,s,1,prec);
    
    arb_clear(T_ratio);
    arb_clear(k_ratio);
    arb_clear(g);
    arb_clear(g_s);
    arb_clear(t);
    arb_clear(s);
    arb_clear(w);
    return 0;
}

void Func_k_to_degrees_of_freedom(arb_t res_r, arb_t res_s, const arb_t k, slong prec)
{
    arb_t kk,T,a,b,error;
    arb_init(kk);
    arb_init(T);
    arb_init(a);
    arb_init(b);
    arb_init(error);
    
    arb_set(kk,k);
    
    //利用 1812.00674 (7) 通过求根的方式来算
    
    arb_set_str(a,"-17",prec);//求根区间，用 ln(T), 可求得 [1，1E18] 内的 k 
    arb_set_str(b,"25",prec);
    arb_set_str(error,"1E-10",prec);
    
    Find_interval_root(T, interior_k_to_degrees_of_freedom, kk, 0,
                       a, b, error,
                       65, Root_Normal, prec);
    
    arb_exp(T,T,prec); // ln(T) --> T
    Effective_degrees_of_freedom_fit(res_r,res_s,T,"Gev",prec);
    
    arb_clear(kk);
    arb_clear(T);
    arb_clear(a);
    arb_clear(b);
    arb_clear(error);
}
