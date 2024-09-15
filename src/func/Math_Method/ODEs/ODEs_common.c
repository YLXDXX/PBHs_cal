#include "../ODEs.h" 
#include <stdlib.h>  

//参考文献 E. Hairer etc.: Solving Ordinary Differential Equations I Nonstiff Problems 「Second Revised Edition」
//计算范数(均方根范数，RMS norm)，用于误差估计和初始步长计算，参看 (4.11)
//注意，在这里我们并没有对每个微分方程分别设置误差，而是统一设定
//进而 sci = Atoli + max(|y0i|, |y1i|) · Rtoli，对每个方程都不同，需传入矢量 sc_i
void ODEs_RMS_norm(arb_t res, const arb_ptr x, const slong dim, const arb_ptr sc_i, slong prec)
{
    arb_t s,sum;
    arb_init(s);
    arb_init(sum);
    
    //arb_dot(s, NULL, 0, x, 1, x, 1, dim, prec);
    arb_zero(sum);
    for(slong i=0; i<dim; i++)
    {
        arb_div(s,x+i,sc_i+i,prec);
        arb_sqr(s,s,prec);
        arb_add(sum,sum,s,prec);
    }
    
    arb_div_si(s,sum,dim,prec);
    arb_sqrt(res,s,prec);
    
    arb_clear(s);
    arb_clear(sum);
}


//自动计算初始步长
void ODEs_select_initial_step(arb_t h, my_odes_func func, const slong dim, void *param, const slong order, //常微分方程组函数
                                const arb_t x_start, const arb_ptr y_start, //给定初始条件
                                const arb_ptr yp_start, //func(x_start,y_start)
                                const arb_t gap_size, //计算区间间隔大小
                                const slong q, //估算方法的阶
                                int direction,
                                const arb_t error_abs, const arb_t error_rel, //绝对误差和相对误差 //误差为绝对精度
                                slong prec)
{
    arb_t s,w,d_0,d_1,d_2,h_0,h_1;
    
    arb_init(s);
    arb_init(w);
    arb_init(d_0);
    arb_init(d_1);
    arb_init(d_2);
    arb_init(h_0);
    arb_init(h_1);
    
    arb_ptr v_s,y_1,sc_i;
    v_s=_arb_vec_init(dim);
    y_1=_arb_vec_init(dim);
    sc_i=_arb_vec_init(dim);
    
    //II.4 Practical Error Estimation and Step Size Selection P169
    //a)
    
    //sc_i== atol + np.abs(y_start) * rtol
    for(slong i=0; i<dim; i++)
    {
        arb_abs(s,y_start+i);
        arb_mul(s,s,error_rel,prec);
        arb_add(sc_i+i,s,error_abs,prec);
    }
    
    
    ODEs_RMS_norm(d_0,y_start,dim,sc_i,prec);
    ODEs_RMS_norm(d_1,yp_start,dim,sc_i,prec);
    
    //b) As a first guess for the step size
    arb_set_si(w,1E5);
    arb_inv(w,w,prec);
    if( arb_lt(d_0,w) || arb_lt(d_1,w) )
    {
        arb_set_si(h_0,1E6);
        arb_inv(h_0,h_0,prec);
    }else
    {
        arb_div(s,d_0,d_1,prec);
        arb_div_si(h_0,s,100,prec);
    }
    
    arb_min(h_0,h_0,gap_size,prec); //不能超过 gap_size
    
    //c) Perform one explicit Euler step
    //y1 = y0 + h0f(x0, y0)
    arb_mul_si(s,h_0,direction,prec); //h0 * direction
    _arb_vec_scalar_mul(y_1,yp_start,dim,s,prec);
    _arb_vec_add(y_1,y_start,y_1,dim,prec);
    
    //f(x0 + h0, y1) 
    arb_add(s,x_start,s,prec); //h0 * direction
    func(v_s, s, y_1, dim, param, order, prec);
    
    //d) Compute d2 
    _arb_vec_sub(v_s,v_s,y_1,dim,prec);
    ODEs_RMS_norm(s,v_s,dim,sc_i,prec);
    arb_div(d_2,s,h_0,prec);
    
    //e) Compute a step size h1
    arb_max(s,d_1,d_2,prec);
    
    arb_one(w);
    arb_div_si(w,w,1E15,prec);
    
    if( arb_lt(s,w) )
    {
        arb_one(s);
        arb_div_si(s,s,1E6,prec);
        arb_one(w);
        arb_div_si(w,w,1E3,prec);
        
        arb_mul(w,h_0,w,prec);
        
        arb_max(h_1,s,w,prec);
    }else
    {
        arb_one(w);
        arb_div_si(w,w,100,prec);
        arb_div(w,w,s,prec);
        
        arb_root_ui(h_1,w,q+1,prec);
    }
    
    //f) Finally we propose as starting step size
    //h = min(100·h0, h1).
    arb_set_si(s,100);
    arb_mul(s,s,h_0,prec);
    
    arb_min(h,s,h_1,prec);
    arb_min(h,h,gap_size,prec); //min(100*h0, h1, gap_size)
    
    arb_clear(s);
    arb_clear(w);
    arb_clear(d_0);
    arb_clear(d_1);
    arb_clear(d_2);
    arb_clear(h_0);
    arb_clear(h_1);
    
    _arb_vec_clear(v_s,dim);
    _arb_vec_clear(y_1,dim);
    _arb_vec_clear(sc_i,dim);
}





