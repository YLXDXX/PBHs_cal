#include "ODEs.h" 
#include <stdlib.h>


/***
int Integration_arb(arb_t res, my_calc_func func, void *param, const slong order,
                    const arb_t a, const arb_t b, const arb_t error,
                    slong step_min , slong step_max,
                    slong prec);

typedef void (*my_odes_func)(arb_ptr yp, const arb_t x, const arb_ptr y, const ulong dim,
                             void* param, const slong order, slong prec);

***/

//求k_1 .. k_6 相关系数
arb_ptr ODEs_RFK45_coe_ah;
arb_ptr ODEs_RFK45_coe_b1,ODEs_RFK45_coe_b2,ODEs_RFK45_coe_b3,ODEs_RFK45_coe_b4,ODEs_RFK45_coe_b5;
arb_ptr ODEs_RFK45_coe_ck_s,ODEs_RFK45_coe_ck_t;

static void Get_RFK45_cal_coe(slong prec)
{
    if(ODEs_RFK45_coe_ah != NULL) //判断是否已设定系数
    {   
        return;
    }
    
    ODEs_RFK45_coe_ah=_arb_vec_init(7);
    ODEs_RFK45_coe_b1=_arb_vec_init(7);
    ODEs_RFK45_coe_b2=_arb_vec_init(7);
    ODEs_RFK45_coe_b3=_arb_vec_init(7);
    ODEs_RFK45_coe_b4=_arb_vec_init(7);
    ODEs_RFK45_coe_b5=_arb_vec_init(7);
    
    ODEs_RFK45_coe_ck_s=_arb_vec_init(7);
    ODEs_RFK45_coe_ck_t=_arb_vec_init(7);
    
    for(int i=0; i<7; i++ )
    {
        arb_one(ODEs_RFK45_coe_ah+i);
        arb_one(ODEs_RFK45_coe_b1+i);
        arb_one(ODEs_RFK45_coe_b2+i);
        arb_one(ODEs_RFK45_coe_b3+i);
        arb_one(ODEs_RFK45_coe_b4+i);
        arb_one(ODEs_RFK45_coe_b5+i);
        arb_one(ODEs_RFK45_coe_ck_s+i);
        arb_one(ODEs_RFK45_coe_ck_t+i);
    }
    
    //ODEs_RFK45_coe_ah k_2 ... k_6
    arb_zero(ODEs_RFK45_coe_ah+0);
    arb_zero(ODEs_RFK45_coe_ah+1);
    
    arb_div_ui(ODEs_RFK45_coe_ah+2,ODEs_RFK45_coe_ah+2,4,prec); //1/4
    
    arb_mul_ui(ODEs_RFK45_coe_ah+3,ODEs_RFK45_coe_ah+3,3,prec); //3/8
    arb_div_ui(ODEs_RFK45_coe_ah+3,ODEs_RFK45_coe_ah+3,8,prec);
    
    arb_mul_ui(ODEs_RFK45_coe_ah+4,ODEs_RFK45_coe_ah+4,12,prec); //12/13
    arb_div_ui(ODEs_RFK45_coe_ah+4,ODEs_RFK45_coe_ah+4,13,prec);
    
    arb_one(ODEs_RFK45_coe_ah+5); //1
    
    arb_div_ui(ODEs_RFK45_coe_ah+6,ODEs_RFK45_coe_ah+6,2,prec); //1/2
    
    
    //ODEs_RFK45_coe_b1 k_2 ... k_6
    arb_zero(ODEs_RFK45_coe_b1+0);
    arb_zero(ODEs_RFK45_coe_b1+1);
    
    arb_div_ui(ODEs_RFK45_coe_b1+2,ODEs_RFK45_coe_b1+2,4,prec); //1/4
    
    arb_mul_ui(ODEs_RFK45_coe_b1+3,ODEs_RFK45_coe_b1+3,3,prec); //3/32
    arb_div_ui(ODEs_RFK45_coe_b1+3,ODEs_RFK45_coe_b1+3,32,prec);
    
    arb_mul_ui(ODEs_RFK45_coe_b1+4,ODEs_RFK45_coe_b1+4,1932,prec); //1932/2197
    arb_div_ui(ODEs_RFK45_coe_b1+4,ODEs_RFK45_coe_b1+4,2197,prec);
    
    arb_mul_ui(ODEs_RFK45_coe_b1+5,ODEs_RFK45_coe_b1+5,439,prec); //439/216
    arb_div_ui(ODEs_RFK45_coe_b1+5,ODEs_RFK45_coe_b1+5,216,prec);
    
    arb_mul_si(ODEs_RFK45_coe_b1+6,ODEs_RFK45_coe_b1+6,-8,prec); //-8/27
    arb_div_ui(ODEs_RFK45_coe_b1+6,ODEs_RFK45_coe_b1+6,27,prec);
    
    
    //ODEs_RFK45_coe_b2 k_3 ... k_6
    arb_zero(ODEs_RFK45_coe_b2+0);
    arb_zero(ODEs_RFK45_coe_b2+1);
    arb_zero(ODEs_RFK45_coe_b2+2);
    
    arb_mul_ui(ODEs_RFK45_coe_b2+3,ODEs_RFK45_coe_b2+3,9,prec); //9/32
    arb_div_ui(ODEs_RFK45_coe_b2+3,ODEs_RFK45_coe_b2+3,32,prec);
    
    arb_mul_si(ODEs_RFK45_coe_b2+4,ODEs_RFK45_coe_b2+4,-7200,prec); //-7200/2197
    arb_div_ui(ODEs_RFK45_coe_b2+4,ODEs_RFK45_coe_b2+4,2197,prec);
    
    arb_mul_si(ODEs_RFK45_coe_b2+5,ODEs_RFK45_coe_b2+5,-8,prec); //-8
    
    arb_mul_si(ODEs_RFK45_coe_b2+6,ODEs_RFK45_coe_b2+6,2,prec); //2
    
    
    //ODEs_RFK45_coe_b3 k_4 ... k_6
    arb_zero(ODEs_RFK45_coe_b3+0);
    arb_zero(ODEs_RFK45_coe_b3+1);
    arb_zero(ODEs_RFK45_coe_b3+2);
    arb_zero(ODEs_RFK45_coe_b3+3);
    
    arb_mul_ui(ODEs_RFK45_coe_b3+4,ODEs_RFK45_coe_b3+4,7296,prec); //7296/2197
    arb_div_ui(ODEs_RFK45_coe_b3+4,ODEs_RFK45_coe_b3+4,2197,prec);
    
    arb_mul_si(ODEs_RFK45_coe_b3+5,ODEs_RFK45_coe_b3+5,3680,prec); //3680/513
    arb_div_ui(ODEs_RFK45_coe_b3+5,ODEs_RFK45_coe_b3+5,513,prec);
    
    arb_mul_si(ODEs_RFK45_coe_b3+6,ODEs_RFK45_coe_b3+6,-3544,prec); //-3544/2565
    arb_div_ui(ODEs_RFK45_coe_b3+6,ODEs_RFK45_coe_b3+6,2565,prec);
    
    //ODEs_RFK45_coe_b4 k_5 ... k_6
    arb_zero(ODEs_RFK45_coe_b4+0);
    arb_zero(ODEs_RFK45_coe_b4+1);
    arb_zero(ODEs_RFK45_coe_b4+2);
    arb_zero(ODEs_RFK45_coe_b4+3);
    arb_zero(ODEs_RFK45_coe_b4+4);
    
    arb_mul_si(ODEs_RFK45_coe_b4+5,ODEs_RFK45_coe_b4+5,-845,prec); //-845/4104
    arb_div_ui(ODEs_RFK45_coe_b4+5,ODEs_RFK45_coe_b4+5,4104,prec);
    
    arb_mul_si(ODEs_RFK45_coe_b4+6,ODEs_RFK45_coe_b4+6,1859,prec); //1859/4104
    arb_div_ui(ODEs_RFK45_coe_b4+6,ODEs_RFK45_coe_b4+6,4104,prec);
    
    
    //ODEs_RFK45_coe_b5 k_6 ... k_6
    arb_zero(ODEs_RFK45_coe_b5+0);
    arb_zero(ODEs_RFK45_coe_b5+1);
    arb_zero(ODEs_RFK45_coe_b5+2);
    arb_zero(ODEs_RFK45_coe_b5+3);
    arb_zero(ODEs_RFK45_coe_b5+4);
    arb_zero(ODEs_RFK45_coe_b5+5);
    
    arb_mul_si(ODEs_RFK45_coe_b5+6,ODEs_RFK45_coe_b5+6,-11,prec); //-11/40
    arb_div_ui(ODEs_RFK45_coe_b5+6,ODEs_RFK45_coe_b5+6,40,prec);
    
    
    //ODEs_RFK45_coe_ck_s  k_1 ... k_5
    arb_zero(ODEs_RFK45_coe_ck_s+0);
    
    arb_mul_si(ODEs_RFK45_coe_ck_s+1,ODEs_RFK45_coe_ck_s+1,25,prec); //25/216
    arb_div_ui(ODEs_RFK45_coe_ck_s+1,ODEs_RFK45_coe_ck_s+1,216,prec);
    
    arb_zero(ODEs_RFK45_coe_ck_s+2); //0
    
    arb_mul_si(ODEs_RFK45_coe_ck_s+3,ODEs_RFK45_coe_ck_s+3,1408,prec); //1408/2565
    arb_div_ui(ODEs_RFK45_coe_ck_s+3,ODEs_RFK45_coe_ck_s+3,2565,prec);
    
    arb_mul_si(ODEs_RFK45_coe_ck_s+4,ODEs_RFK45_coe_ck_s+4,2197,prec); //2197/4104
    arb_div_ui(ODEs_RFK45_coe_ck_s+4,ODEs_RFK45_coe_ck_s+4,4104,prec);
    
    arb_mul_si(ODEs_RFK45_coe_ck_s+5,ODEs_RFK45_coe_ck_s+5,-1,prec); //-1/5
    arb_div_ui(ODEs_RFK45_coe_ck_s+5,ODEs_RFK45_coe_ck_s+5,5,prec);
    
    arb_zero(ODEs_RFK45_coe_ck_s+6); //0
    
    
    //ODEs_RFK45_coe_ck_t  k_1 ... k_5
    arb_zero(ODEs_RFK45_coe_ck_t+0);
    
    arb_mul_si(ODEs_RFK45_coe_ck_t+1,ODEs_RFK45_coe_ck_t+1,16,prec); //16/135
    arb_div_ui(ODEs_RFK45_coe_ck_t+1,ODEs_RFK45_coe_ck_t+1,135,prec);
    
    arb_zero(ODEs_RFK45_coe_ck_t+2); //0
    
    arb_mul_si(ODEs_RFK45_coe_ck_t+3,ODEs_RFK45_coe_ck_t+3,6656,prec); //6656/12825
    arb_div_ui(ODEs_RFK45_coe_ck_t+3,ODEs_RFK45_coe_ck_t+3,12825,prec);
    
    arb_mul_si(ODEs_RFK45_coe_ck_t+4,ODEs_RFK45_coe_ck_t+4,28561,prec); //28561/56430
    arb_div_ui(ODEs_RFK45_coe_ck_t+4,ODEs_RFK45_coe_ck_t+4,56430,prec);
    
    arb_mul_si(ODEs_RFK45_coe_ck_t+5,ODEs_RFK45_coe_ck_t+5,-9,prec); //-9/50
    arb_div_ui(ODEs_RFK45_coe_ck_t+5,ODEs_RFK45_coe_ck_t+5,50,prec);
    
    arb_mul_si(ODEs_RFK45_coe_ck_t+6,ODEs_RFK45_coe_ck_t+6,2,prec); //2/55
    arb_div_ui(ODEs_RFK45_coe_ck_t+6,ODEs_RFK45_coe_ck_t+6,55,prec);
    
}

int ODEs_RFK45(arb_ptr res, my_odes_func func, const slong dim, void *param, const slong order, //常微分方程组函数
                    const arb_t x_start, const arb_ptr y_start, //给定初始条件
                    const arb_t x_end, //求出点 x_end 对应的函数值
                    const slong num, const arb_t error, //num为迭代区间为[x_start,x_end]，将其分为几等份，从而给出初始步长
                                                        //误差为相对精度
                    slong prec)
{
    Get_RFK45_cal_coe(prec);
    
    arb_t s,t,h,x,x_tem,pow_1_div_4,coe_0_dot_84,delta;
    arb_init(s);
    arb_init(t);
    arb_init(h);
    arb_init(x);
    arb_init(x_tem);
    arb_init(pow_1_div_4);
    arb_init(coe_0_dot_84);
    arb_init(delta);
    
    arb_one(pow_1_div_4); //1/4
    arb_div_ui(pow_1_div_4,pow_1_div_4,4,prec);
    
    arb_one(coe_0_dot_84); //0.84
    arb_mul_ui(coe_0_dot_84,coe_0_dot_84,84,prec);
    arb_div_ui(coe_0_dot_84,coe_0_dot_84,100,prec);
    
    
    arb_ptr w,y_i,y_i1_tt,y_i1_ss,y_tem,k1,k2,k3,k4,k5,k6,R;
    w=_arb_vec_init(dim);
    y_i=_arb_vec_init(dim);
    y_i1_tt=_arb_vec_init(dim);
    y_i1_ss=_arb_vec_init(dim);
    y_tem=_arb_vec_init(dim);
    
    k1=_arb_vec_init(dim);
    k2=_arb_vec_init(dim);
    k3=_arb_vec_init(dim);
    k4=_arb_vec_init(dim);
    k5=_arb_vec_init(dim);
    k6=_arb_vec_init(dim);
    R=_arb_vec_init(dim);
    
    
    _arb_vec_set(y_i,y_start,dim); //初始条件
    
    
    //迭代初始步长
    arb_sub(h,x_end,x_start,prec);
    arb_div_ui(h,h,num,prec);
    
    
    int judge_toward,judge_finish;
    
    judge_toward=0;
    judge_finish=1;
    
    if( arb_is_negative(h) )
    {
        judge_toward=1;
    }
        
    arb_set(x,x_start); //初始迭代位置
    
    for(slong i=0; ; i++) //这里会自动调节步长，不能使用 i<num 来判断退出迭代
    {
        //相关系数求解
        //k1
        func(k1, x, y_i, dim, param, order, prec);
        _arb_vec_scalar_mul(k1,k1,dim,h,prec);
        
        //k2
        arb_mul(x_tem,ODEs_RFK45_coe_ah+2,h,prec);
        arb_add(x_tem,x_tem,x,prec);
        
        _arb_vec_scalar_mul(y_tem,k1,dim,ODEs_RFK45_coe_b1+2,prec);
        _arb_vec_add(y_tem,y_tem,y_i,dim,prec);
        
        func(k2, x_tem, y_tem, dim, param, order, prec);
        _arb_vec_scalar_mul(k2,k2,dim,h,prec);
        
        //k3
        arb_mul(x_tem,ODEs_RFK45_coe_ah+3,h,prec);
        arb_add(x_tem,x_tem,x,prec);
        
        _arb_vec_scalar_mul(y_tem,k1,dim,ODEs_RFK45_coe_b1+3,prec);
        _arb_vec_add(y_tem,y_tem,y_i,dim,prec);
        _arb_vec_scalar_mul(w,k2,dim,ODEs_RFK45_coe_b2+3,prec);
        _arb_vec_add(y_tem,y_tem,w,dim,prec);
        
        func(k3, x_tem, y_tem, dim, param, order, prec);
        _arb_vec_scalar_mul(k3,k3,dim,h,prec);
        
        
        //k4
        arb_mul(x_tem,ODEs_RFK45_coe_ah+4,h,prec);
        arb_add(x_tem,x_tem,x,prec);
        
        _arb_vec_scalar_mul(y_tem,k1,dim,ODEs_RFK45_coe_b1+4,prec);
        _arb_vec_add(y_tem,y_tem,y_i,dim,prec);
        _arb_vec_scalar_mul(w,k2,dim,ODEs_RFK45_coe_b2+4,prec);
        _arb_vec_add(y_tem,y_tem,w,dim,prec);
        _arb_vec_scalar_mul(w,k3,dim,ODEs_RFK45_coe_b3+4,prec);
        _arb_vec_add(y_tem,y_tem,w,dim,prec);
        
        func(k4, x_tem, y_tem, dim, param, order, prec);
        _arb_vec_scalar_mul(k4,k4,dim,h,prec);
        
        //k5
        arb_mul(x_tem,ODEs_RFK45_coe_ah+5,h,prec);
        arb_add(x_tem,x_tem,x,prec);
        
        _arb_vec_scalar_mul(y_tem,k1,dim,ODEs_RFK45_coe_b1+5,prec);
        _arb_vec_add(y_tem,y_tem,y_i,dim,prec);
        _arb_vec_scalar_mul(w,k2,dim,ODEs_RFK45_coe_b2+5,prec);
        _arb_vec_add(y_tem,y_tem,w,dim,prec);
        _arb_vec_scalar_mul(w,k3,dim,ODEs_RFK45_coe_b3+5,prec);
        _arb_vec_add(y_tem,y_tem,w,dim,prec);
        _arb_vec_scalar_mul(w,k4,dim,ODEs_RFK45_coe_b4+5,prec);
        _arb_vec_add(y_tem,y_tem,w,dim,prec);
        
        func(k5, x_tem, y_tem, dim, param, order, prec);
        _arb_vec_scalar_mul(k5,k5,dim,h,prec);
        
        
        //k6
        arb_mul(x_tem,ODEs_RFK45_coe_ah+6,h,prec);
        arb_add(x_tem,x_tem,x,prec);
        
        _arb_vec_scalar_mul(y_tem,k1,dim,ODEs_RFK45_coe_b1+6,prec);
        _arb_vec_add(y_tem,y_tem,y_i,dim,prec);
        _arb_vec_scalar_mul(w,k2,dim,ODEs_RFK45_coe_b2+6,prec);
        _arb_vec_add(y_tem,y_tem,w,dim,prec);
        _arb_vec_scalar_mul(w,k3,dim,ODEs_RFK45_coe_b3+6,prec);
        _arb_vec_add(y_tem,y_tem,w,dim,prec);
        _arb_vec_scalar_mul(w,k4,dim,ODEs_RFK45_coe_b4+6,prec);
        _arb_vec_add(y_tem,y_tem,w,dim,prec);
        _arb_vec_scalar_mul(w,k5,dim,ODEs_RFK45_coe_b5+6,prec);
        _arb_vec_add(y_tem,y_tem,w,dim,prec);
        
        func(k6, x_tem, y_tem, dim, param, order, prec);
        _arb_vec_scalar_mul(k6,k6,dim,h,prec);
        
        //对于y_{i+1}的两个估计值
        _arb_vec_scalar_mul(y_tem,k1,dim,ODEs_RFK45_coe_ck_s+1,prec);
        _arb_vec_add(y_tem,y_tem,y_i,dim,prec);
        _arb_vec_scalar_mul(w,k3,dim,ODEs_RFK45_coe_ck_s+3,prec);
        _arb_vec_add(y_tem,y_tem,w,dim,prec);
        _arb_vec_scalar_mul(w,k4,dim,ODEs_RFK45_coe_ck_s+4,prec);
        _arb_vec_add(y_tem,y_tem,w,dim,prec);
        _arb_vec_scalar_mul(w,k5,dim,ODEs_RFK45_coe_ck_s+5,prec);
        _arb_vec_add(y_i1_ss,y_tem,w,dim,prec);
        
        
        _arb_vec_scalar_mul(y_tem,k1,dim,ODEs_RFK45_coe_ck_t+1,prec);
        _arb_vec_add(y_tem,y_tem,y_i,dim,prec);
        _arb_vec_scalar_mul(w,k3,dim,ODEs_RFK45_coe_ck_t+3,prec);
        _arb_vec_add(y_tem,y_tem,w,dim,prec);
        _arb_vec_scalar_mul(w,k4,dim,ODEs_RFK45_coe_ck_t+4,prec);
        _arb_vec_add(y_tem,y_tem,w,dim,prec);
        _arb_vec_scalar_mul(w,k5,dim,ODEs_RFK45_coe_ck_t+5,prec);
        _arb_vec_add(y_tem,y_tem,w,dim,prec);
        _arb_vec_scalar_mul(w,k6,dim,ODEs_RFK45_coe_ck_t+6,prec);
        _arb_vec_add(y_i1_tt,y_tem,w,dim,prec);
        
        
        //得到计算的误差
        _arb_vec_sub(R,y_i1_ss,y_i1_tt,dim,prec);
        _arb_vec_scalar_div(R,R,dim,h,prec);
        
        
        for(slong j=0; j< dim;j++) //在R中筛选出最大的元素，用以后面误差估算 
        {
            arb_abs(R+j,R+j); //每组的计算误差
            
            if( arb_gt(R+j,R) )
            {
                arb_swap(R+j,R); //这里，保证 R 中第一个元素是所有元素中最大的那个
            }
        }
        
        arb_div(s,error,R,prec);
        arb_pow(s,s,pow_1_div_4,prec);
        arb_mul(delta,s,coe_0_dot_84,prec);
        arb_get_mid_arb(delta,delta); //消除每次迭代的误差累积
        
        if( arb_lt(R,error) ) //到精度要求
        {
            if(judge_finish==0)
            {
                //完成相应的计算
                _arb_vec_set(res,y_i1_ss,dim);
                break;
            }
            
            arb_add(x,x,h,prec); //改变迭代位置，这里需使用旧步长
            _arb_vec_swap(y_i1_ss,y_i,dim); //同时改变迭代初始条件
        }
        
        
        arb_mul(h,delta,h,prec); //改变下一次迭代的步长，需在X位置更新后
        
        arb_sub(t,x_end,x,prec); //用于终点
        
        if( judge_toward==0 && arb_gt(h,t))
        {
            arb_set(h,t);
            judge_finish=0;
        }
        
        if(judge_toward==1 && arb_lt(h,t))
        {
            arb_set(h,t);
            judge_finish=0;
        }
        //printf("%li\n\n",i);
    }
    
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(h);
    arb_clear(x);
    arb_clear(x_tem);
    arb_clear(pow_1_div_4);
    arb_clear(coe_0_dot_84);
    arb_clear(delta);
    
    _arb_vec_clear(w,dim);
    _arb_vec_clear(y_i,dim);
    _arb_vec_clear(y_i1_tt,dim);
    _arb_vec_clear(y_i1_ss,dim);
    _arb_vec_clear(y_tem,dim);
    
    _arb_vec_clear(k1,dim);
    _arb_vec_clear(k2,dim);
    _arb_vec_clear(k3,dim);
    _arb_vec_clear(k4,dim);
    _arb_vec_clear(k5,dim);
    _arb_vec_clear(k6,dim);
    _arb_vec_clear(R,dim);
    
    
    return 0;
}


