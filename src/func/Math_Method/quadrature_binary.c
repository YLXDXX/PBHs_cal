#include "quadrature_binary.h"

//二元积分，积分区域为矩形
int integration_binary_rectangle(arb_t res, my_calc_func_binary func, void *param, const slong order,
                                 const arb_t x_a, const arb_t x_b, const arb_t x_error, 
                                 slong x_step_min , slong x_step_max,
                                 const arb_t y_a, const arb_t y_b, const arb_t y_error,
                                 slong y_step_min , slong y_step_max,
                                 slong prec)
{
    arb_t t;
    arb_init(t);
    
    int res_judge=0;
    
    //根据需要采用不同的积分方法
    switch(Integral_method)
    {
        case gauss_kronrod_iterate : //gauss_kronrod迭代版本
            get_gauss_kronrod_node_weight(65,prec); //获取 gauss_kronrod 的节点位置和权重
            res_judge=integration_binary_rectangle_gauss_kronrod(t, func, param, order,
                                                                 x_a, x_b, x_error,
                                                                 x_step_min, x_step_max,
                                                                 y_a, y_b, y_error,
                                                                 y_step_min, y_step_max,
                                                                 prec);
            break;
            
        case double_exponential :
            res_judge=quadrature_binary_rectangle_double_exponential(t, func, param, order,
                                                                     x_a, x_b, x_error,
                                                                     x_step_min, x_step_max,
                                                                     y_a, y_b, y_error,
                                                                     y_step_min, y_step_max,
                                                                     prec);
            break;
            
        default : //默认使用gauss_kronrod_iterate积分
            get_gauss_kronrod_node_weight(65,prec); //获取 gauss_kronrod 的节点位置和权重
            res_judge=integration_binary_rectangle_gauss_kronrod(t, func, param, order,
                                                                 x_a, x_b, x_error,
                                                                 x_step_min, x_step_max,
                                                                 y_a, y_b, y_error,
                                                                 y_step_min, y_step_max,
                                                                 prec);
            
    }
    
    arb_set(res,t);
    arb_clear(t);
    
    if(res_judge==0)
    {
        return 0;
    }else
    {
        return 1;
    }
    
}


//二元积分，积分区域为矩形，对矩形区域自适应
int integration_binary_rectangle_adaptive(arb_t res, my_calc_func_binary func, void *param, const slong order,
                                          const arb_t x_a, const arb_t x_b, const arb_t x_error, 
                                          slong x_step_min , slong x_step_max,
                                          const arb_t y_a, const arb_t y_b, const arb_t y_error,
                                          slong y_step_min , slong y_step_max,
                                          slong prec)
{
    arb_t t;
    arb_init(t);
    
    int res_judge=0;
    
    //根据需要采用不同的积分方法
    switch(Integral_method)
    {
        case gauss_kronrod_iterate : //gauss_kronrod迭代版本
            get_gauss_kronrod_node_weight(65,prec); //获取 gauss_kronrod 的节点位置和权重
            res_judge=integration_binary_rectangle_adaptive_gauss_kronrod(t, func, param, order,
                                                                 x_a, x_b, x_error,
                                                                 x_step_min, x_step_max,
                                                                 y_a, y_b, y_error,
                                                                 y_step_min, y_step_max,
                                                                 prec);
            break;
            
        case double_exponential :
            res_judge=integration_binary_rectangle_adaptive_double_exponential(t, func, param, order,
                                                                     x_a, x_b, x_error,
                                                                     x_step_min, x_step_max,
                                                                     y_a, y_b, y_error,
                                                                     y_step_min, y_step_max,
                                                                     prec);
            break;
            
        default : //默认使用gauss_kronrod_iterate积分
            get_gauss_kronrod_node_weight(65,prec); //获取 gauss_kronrod 的节点位置和权重
            res_judge=integration_binary_rectangle_adaptive_gauss_kronrod(t, func, param, order,
                                                                 x_a, x_b, x_error,
                                                                 x_step_min, x_step_max,
                                                                 y_a, y_b, y_error,
                                                                 y_step_min, y_step_max,
                                                                 prec);
            
    }
    
    arb_set(res,t);
    arb_clear(t);
    
    if(res_judge==0)
    {
        return 0;
    }else
    {
        return 1;
    } 
}

//二元积分，积分区域非矩形
int integration_binary_func(arb_t res, my_calc_func_binary func, void* param, const slong order,
                            const arb_t x_a, const arb_t x_b, const arb_t x_error, 
                            slong x_step_min , slong x_step_max,
                            my_calc_func y_a_func, void *param_y_a,  const slong order_y_a,
                            my_calc_func y_b_func, void *param_y_b,  const slong order_y_b,
                            const arb_t y_error, slong y_step_min, slong y_step_max,
                            slong prec)
{
    arb_t t;
    arb_init(t);
    
    int res_judge=0;
    
    //根据需要采用不同的积分方法
    switch(Integral_method)
    {
        case gauss_kronrod_iterate : //gauss_kronrod迭代版本
            get_gauss_kronrod_node_weight(65,prec); //获取 gauss_kronrod 的节点位置和权重
            res_judge=integration_binary_func_gauss_kronrod(t, func, param, order,
                                                            x_a, x_b, x_error,
                                                            x_step_min, x_step_max,
                                                            y_a_func, param_y_a, order_y_a,
                                                            y_b_func, param_y_b, order_y_b,
                                                            y_error, y_step_min, y_step_max,
                                                            prec);
            break;
            
        case double_exponential :
            res_judge=integration_binary_func_double_exponential(t, func, param, order,
                                                                 x_a, x_b, x_error,
                                                                 x_step_min, x_step_max,
                                                                 y_a_func, param_y_a, order_y_a,
                                                                 y_b_func, param_y_b, order_y_b,
                                                                 y_error, y_step_min, y_step_max,
                                                                 prec);
            break;
            
        default : //默认使用gauss_kronrod_iterate积分
            get_gauss_kronrod_node_weight(65,prec); //获取 gauss_kronrod 的节点位置和权重
            res_judge=integration_binary_func_gauss_kronrod(t, func, param, order,
                                                            x_a, x_b, x_error,
                                                            x_step_min, x_step_max,
                                                            y_a_func, param_y_a, order_y_a,
                                                            y_b_func, param_y_b, order_y_b,
                                                            y_error, y_step_min, y_step_max,
                                                            prec);
            
    }
    
    arb_set(res,t);
    arb_clear(t);
    
    if(res_judge==0)
    {
        return 0;
    }else
    {
        return 1;
    }
}
