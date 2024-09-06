#include "routine_test.h"
#include <stdlib.h>

//测试程序
void routine_test(slong prec)
{
    //微分方程测试
    
    arb_t s,t,w,x_start,x_end,error;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(x_start);
    arb_init(x_end);
    arb_init(error);
    
    int dim=2;
    
    arb_ptr v_s,v_t,y_start;
    v_s=_arb_vec_init(dim);
    v_t=_arb_vec_init(dim);
    y_start=_arb_vec_init(dim);
    
    arb_set_str(x_start,"0",prec); //初始条件
    arb_set_str(y_start,"1",prec);
    arb_set_str(y_start+1,"1",prec);
    arb_set_str(error,"1E-8",prec);
    
    arb_set_str(x_end,"-15",prec);
    
    
    
    ODEs_RFK45(v_s, Func_test_ODEs_func_04, dim, NULL, 0, //常微分方程组函数
               x_start, y_start, //给定初始条件
               x_end, //求出点 x_end 对应的函数值
               50, error, //num为迭代区间为[x_start,x_end]，将其分为几等份，从而给出初始步长
               prec);
    
    arb_printn(v_s, 50,0);printf("\n");
    
    
    //arb_set_str(x_end,"0.022",prec);
    
    /***
    // Func_test_ODEs_func_01 精确解为 y=1/x^3 * Exp(1-1/x)
    arb_pow_ui(s,x_end,3,prec);
    arb_inv(t,x_end,prec);
    arb_neg(t,t);
    arb_add_ui(t,t,1,prec);
    arb_exp(t,t,prec);
    arb_div(s,t,s,prec);
    ***/
    
    /***
    // Func_test_ODEs_func_02 精确解为 y=x^2+2*x+1-e^x/2
    arb_sqr(s,x_end,prec);
    arb_mul_si(t,x_end,2,prec);
    arb_add(s,s,t,prec);
    arb_add_ui(s,s,1,prec);
    arb_exp(t,x_end,prec);
    arb_div_ui(t,t,2,prec);
    arb_sub(s,s,t,prec);
    ***/
    
    /*
    // Func_test_ODEs_func_03 精确解为 y=cos(x)+sin(x)+x^2-2
    arb_sin(s,x_end,prec);
    arb_cos(t,x_end,prec);
    arb_add(s,s,t,prec);
    arb_sqr(t,x_end,prec);
    arb_add(s,s,t,prec);
    arb_sub_ui(s,s,2,prec);
    */
    
    // Func_test_ODEs_func_04 精确解为 y=(1+1/6*x^3)*e^{x}
    
    arb_pow_ui(s,x_end,3,prec);
    arb_div_ui(s,s,6,prec);
    arb_add_ui(s,s,1,prec);
    
    arb_exp(t,x_end,prec);
    arb_mul(s,s,t,prec);
    
    arb_printn(s, 50,0);printf("\n");
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    
    _arb_vec_clear(v_s,dim);
    _arb_vec_clear(v_t,dim);
    
    exit(0);
    
    //积分测试 1
    
    arb_t x,y,x_a,x_b,y_a,y_b,e,r,k,eta;
    
    arb_init(x);
    arb_init(y);
    arb_init(x_a);
    arb_init(x_b);
    arb_init(y_a);
    arb_init(y_b);
    arb_init(e);
    arb_init(r);
    arb_init(k);
    arb_init(eta);
    
    Integral_method=gauss_kronrod_iterate; // gauss_kronrod_iterate/double_exponential
    
    //诱导引力波积分辐助函数的积分设定
    arb_zero(Int_GW_I_func_min); //积分下限，为 1/k or 0 or 1，这里我们设为1
    arb_pos_inf(Int_GW_I_func_max); //积分上限 +∞
    //arb_set_str(Int_GW_I_func_max,"1E4",prec);
    arb_set_str(Int_GW_I_func_precision,"1E-10",prec);
    Int_GW_I_func_iterate_min=6;
    Int_GW_I_func_iterate_max=15;
        
    
    //诱导引力波功率谱的积分设定
    arb_set_str(Int_GW_power_spectra_min,"1E-15",prec);
    arb_set_str(Int_GW_power_spectra_max,"1E3",prec);
    arb_set_str(Int_GW_power_spectra_precision,"1E-15",prec);
    Int_GW_power_spectra_iterate_min=6;
    Int_GW_power_spectra_iterate_max=15;
    
    // 功率谱积矩形分用
    arb_set_str(Int_GW_power_spectra_x_min,"0",prec); //对x积分
    arb_set_str(Int_GW_power_spectra_x_max,"1E3",prec);
    arb_set_str(Int_GW_power_spectra_x_precision,"1E-15",prec);
    Int_GW_power_spectra_iterate_x_min=4;
    Int_GW_power_spectra_iterate_x_max=15;
    
    arb_set_str(Int_GW_power_spectra_y_min,"-1",prec); //对y积分
    arb_set_str(Int_GW_power_spectra_y_max,"1",prec);
    arb_set_str(Int_GW_power_spectra_y_precision,"1E-15",prec);
    Int_GW_power_spectra_iterate_y_min=6;
    Int_GW_power_spectra_iterate_y_max=100;
    
    /*
    arb_set_str(x,"9E6",prec);
    arb_log(x,x,prec);
    power_spectrum(y,x,prec);
    
    arb_printn(y, 50,0);printf("\n");
    exit(0);
    */
    
    /*
    //诱导引力波
    arb_set_str(x,"2",prec);
    arb_set_str(y,"2.1",prec);
    
    GW_I_s_func_analyze(r,x,y,x,prec);
    arb_printn(r,50,0);printf("\n");
    
    
    GW_I_s_func(r,x,y,x,prec);
    arb_printn(r,50,0);printf("\n");
    
    
    Integral_method=gauss_kronrod_iterate;Int_GW_I_func_iterate_min=128;Int_GW_I_func_iterate_max=10000;
    GW_I_s_func(r,x,y,x,prec);
    arb_printn(r,50,0);printf("\n");
    
    exit(0);
    */
    
    
    arb_set_str(k,"1.56E6",prec);
    //arb_set_str(eta,"1",prec);
    arb_inv(eta,k,prec);
    //GW_power_spectra(r,eta,k,prec);
    //GW_current_energy_density_01(r,k,prec);
    //GW_current_energy_density_02(r,k,prec);
    //GW_current_energy_density_Omega_G(r,k,prec);
    //GW_current_energy_density_Omega_dim_8(r,x,k,prec);
    //GW_current_energy_density_Omega_dim_6(r,k,prec);
    //GW_current_energy_density_Omega_dim_5(r,x,k,prec);
    //GW_current_energy_density_Omega_dim_4(r,k,prec);
    //GW_current_energy_density_Omega_dim_2(r,k,prec);
    //GW_current_energy_density(r,k,prec);
    //GW_current_energy_density_cuba(r,k,0,prec);
    arb_printn(r,50,0);printf("\n");
    arb_printn(x,50,0);printf("\n");
    
    //画图
    //draw_pic(NULL,prec); //输出点用于画图，可从命令行传递参数
    
    exit(0);
    //GW_power_spectra_2(x,eta,k,prec);
    GW_current_energy_density_02(x,k,prec);
    arb_div(r,x,r,prec);
    arb_printn(r,50,0);printf("\n");
    exit(0);
    
    
    //验证概率守恒
    /*
    arb_t s,a,b;
    arb_init(s);
    arb_init(a);
    arb_init(b);
    
    arb_set_str(a,"-1.5",prec); //这里 -infinity 取1
    
    arb_abs(b,Up_step_h);
    arb_inv(b,b,prec);
    arb_mul_si(b,b,2,prec);
    
    arb_set_str(s,"1E-30",prec); //精度
    integration_gauss_kronrod(s, Probability_zeta, NULL, 0, 
                              a, b,s,
                              Integration_iterate_min,Integration_iterate_max, prec);
    arb_printn(s,25,0);printf("\n");
    */
    
    
    /*
    //积分测试 1
    arb_set_str(x_a,"0.5E20",prec);
    arb_set_str(x_b,"5E20",prec);
    arb_set_str(y_a,"0.5E20",prec);
    arb_set_str(y_b,"5E20",prec);
    arb_set_str(e,"1E-10",prec);
    integration_binary(r,Func_test_3,NULL,0,
                       x_a,x_b,e,5,500,
                       y_a,y_b,e,5,500,
                       prec);
    arb_printn(r, 30, 0);printf("\n");
    */
    
    /*
    //积分测试 2
    arb_set_str(x_a,"0",prec);
    arb_set_str(x_b,"18",prec);
    arb_set_str(e,"1E-10",prec);
    Integration_arb(r, Func_test, NULL, 0, 
                              x_a, x_b,e,
                              Integration_iterate_min,Integration_iterate_max, prec);
    arb_printn(r, 50,0);printf("\n");
    */
    
    
    /*
    //积分测试 3
    arb_set_str(x_a,"0",prec);
    arb_set_str(x_b,"498",prec);
    //arb_neg_inf(x_a);
    //arb_pos_inf(x_b);
    arb_set_str(e,"1E-30",prec);
    Double_Exponential_Quadrature(r,Func_test_9,NULL,0,
                                  x_a,x_b,e,5,16,prec);
    arb_printn(r, 50, 0);printf("\n");
    Integration_arb(r, Func_test_9, NULL, 0, 
                              x_a, x_b,e,
                              Integration_iterate_min,Integration_iterate_max, prec);
    arb_printn(r, 50, 0);printf("\n");
    exit(0);
    */
    
    //发散函数积分
    arb_set_str(x_a,"-10",prec);
    arb_set_str(x_b,"10",prec);
    arb_set_str(y_a,"30",prec);
    arb_set_str(y_b,"-30",prec);
    
    //积分行为设定
    //Integral_method=gauss_kronrod_iterate; // gauss_kronrod_iterate/double_exponential
    
    //arb_neg_inf(x_a);
    //arb_pos_inf(x_b);
    //arb_zero(x_b);
    
    arb_set_str(e,"1E-10",prec);
    
    /*
    //测试二重积分，积分区域为非矩形
    integration_binary_func(r, Func_test_quad_func_03, NULL, 0,
                            x_a, x_b, e, 
                            2, 13,
                            Func_test_quad_func_03_y_a, NULL,  0,
                            Func_test_quad_func_03_y_b, NULL,  0,
                            e, 2, 12,
                            prec);
    */
    
    //测试二重积分，积分区域为矩形
    
    Integral_method=double_exponential;
    integration_binary_rectangle_adaptive(r, Func_test_quad_rectangle_01, NULL, 0,
                                          x_a, x_b, e, 
                                          2, 18,
                                          y_a, y_b, e, 2, 18,
                                          prec);
    arb_printn(r, 50, 0);printf("\n");
    Integral_method=gauss_kronrod_iterate;
    integration_binary_rectangle_adaptive(r, Func_test_quad_rectangle_01, NULL, 0,
                                          x_a, x_b, e, 
                                          4, 1300,
                                          y_a, y_b, e, 4, 1300,
                                          prec);
    arb_printn(r, 50, 0);printf("\n");
    exit(0);
    
    
    //arb_neg_inf(x_a);
    //arb_pos_inf(x_b);
    arb_set_str(e,"1E-30",prec);
    //注意到，这里的结果应为无穷，其积分值受精度影响，精度越高，积分值越大
    Double_Exponential_Quadrature(r,Func_test_8,NULL,0,
                                  x_a,x_b,e,5,16,prec);
    arb_printn(r, 50, 0);printf("\n");
    Integration_arb(r, Func_test_8, NULL, 0, 
                    x_a, x_b,e,5,160, prec);
    arb_printn(r, 50, 0);printf("\n");
    exit(0);
    
    
    /*
    //找根测试
    arb_set_str(x_a,"0.01",prec);
    arb_set_str(x_b,"80",prec);
    arb_set_str(e,"1E-10",prec);
    Find_interval_root(r,Find_root_test,NULL,0,x_a,x_b,e,10,prec);
    arb_printn(r, 30, 0);printf("\n");
    */
    
    
    
    /*
    //多根测试
    arb_set_str(x_a,"-2",prec);
    arb_set_str(x_b,"+10",prec);
    arb_set_str(e,"1E-13",prec);
    arb_ptr muil_r; //muil_r是一个指针，
    arb_ptr* m_r; //改变muil_r指针指向的地址，需要一个指向该指针的指针
    m_r=&muil_r;
    //muil_r=_arb_vec_init(2);
    int root_num;
    root_num=Find_interval_multi_root(m_r,Func_test_4,NULL, 0,
                             x_a,x_b,e,6500,prec);
    printf("根的个数为: %i\n",root_num);
    for(int r_i=0; r_i<root_num; r_i++)
    {
        arb_printn(muil_r+r_i, 50,0);printf("\n");
    }
    _arb_vec_clear(muil_r, root_num); //清理数组
    */
    
    arb_clear(y);
    arb_clear(x);
    arb_clear(x_a);
    arb_clear(x_b);
    arb_clear(y_a);
    arb_clear(y_b);
    arb_clear(e);
    arb_clear(r);
}
 
