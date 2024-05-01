#include "routine_test.h"
#include <stdlib.h>

//测试程序
void routine_test(slong prec)
{
    //积分测试 1
    
    arb_t x_a,x_b,y_a,y_b,e,r;
    
    arb_init(x_a);
    arb_init(x_b);
    arb_init(y_a);
    arb_init(y_b);
    arb_init(e);
    arb_init(r);
    
    
    
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
    arb_set_str(x_a,"0",prec);
    arb_set_str(x_b,"1",prec);
    //arb_neg_inf(x_a);
    //arb_pos_inf(x_b);
    arb_set_str(e,"1E-100",prec);
    //注意到，这里的结果应为无穷，其积分值受精度影响，精度越高，积分值越大
    Double_Exponential_Quadrature(r,Func_test_10,NULL,0,
                                  x_a,x_b,e,5,16,prec);
    arb_printn(r, 50, 0);printf("\n");
    Integration_arb(r, Func_test_10, NULL, 0, 
                    x_a, x_b,e,
                    32,10000, prec);
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
    
    arb_clear(x_a);
    arb_clear(x_b);
    arb_clear(y_a);
    arb_clear(y_b);
    arb_clear(e);
    arb_clear(r);
}
 
