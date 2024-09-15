#include "../ODEs.h" 
#include <stdlib.h>


//Runge–Kutta 中的 RFK45 方法，嵌入方法 4(5), 阶数 q=4，用于误差估算阶数 q=5 
//暂无 dense output

//为了得到 RFK45 求解中间各点的值，专门定义了一个结构体来处理
//主要是利用插值来处理，但是该方法有问题，当求解跨度特别大，如 x_end-x_start=1E6 时
//并不能跟随自动迭代步骤，其计算点一多，耗时很长，精度损失严重
//后续可使用 dense output 方法优化（如 Hermite 三次多项式插值）
//这里的点输出，需要事先指定输出那些点

//初始化微分方程多点输出结构，为其分配内存
ODEs_point_output_t ODEs_point_output_init(slong num, slong dim)
{
    //通过 ODEs_point_output_t 声明变量后，需将该结构体指针指向结构体
    ODEs_point_output_t p_out= (ODEs_point_output_t)calloc(1,sizeof(struct ODEs_point_output_structure));
    
    
    p_out->num=num;
    p_out->dim=dim;
    
    p_out->p_x=_arb_vec_init(num);
    p_out->p_y=(arb_ptr*)calloc(num,sizeof(arb_ptr));
    
    for(slong i=0; i<num; i++) //对每个系数的多项式表达进行初始化
    {
        p_out->p_y[i]=_arb_vec_init(dim);
    }
    
    return p_out;
}


//清理微分方程多点输出结构，释放内存
void ODEs_point_output_clear(ODEs_point_output_t p_out, slong num, slong dim)
{
    _arb_vec_clear(p_out->p_x,num);
    
    //首先释放每个多项式系数
    for(slong i=0; i<num; i++)
    {
        _arb_vec_clear(p_out->p_y[i],dim);
    }
    
    free(p_out);
}


// RFK45 计算主程序
//这里，没有使用区分相对误差和绝对误差，没有使用自动步长，没有 dense output
//后续可优化，参照 DOPRI54
int ODEs_RFK45(arb_ptr y_end, my_odes_func func, const slong dim, void *param, const slong order, //常微分方程组函数
                    const arb_t x_start, const arb_ptr y_start, //给定初始条件
                    const arb_t x_end, //求出点 x_end 对应的函数值
                    const slong num, const arb_t error, //num为迭代区间为[x_start,x_end]，将其分为几等份，从而给出初始步长
                                                        //误差为相对精度
                    slong prec)
{
    //Get_RFK45_cal_coe(prec); //获取系数
    //计算相关系数为全局变量，在程序运行前，在全局变量部分进行初始化
    
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
    if( num<2 ) //注意，这里 num 不能为 1，否则在会面取 h 时，会刚好抵到 x_end , 导致错误
    {
        arb_sub(h,x_end,x_start,prec);
        arb_div_ui(h,h,2,prec);
    }else
    {
        arb_sub(h,x_end,x_start,prec);
        arb_div_ui(h,h,num,prec);
    }
    
    
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
            
            arb_get_mid_arb(y_i1_ss+j,y_i1_ss+j); //迭代次数足够多时，计算误差会很大
                                                //为了不影响后面计算，如 ODEs_RFK45_interval_point_output，这里进行了消除
        }
        
        if(!arb_is_finite(R))
        {
            printf("遇到错误已退出，可能原因：\n");
            printf("1) 精度要求太高，当前精度不能满足要求，请提高 prec 值后重试\n");
            printf("2) 求解条件设定有问题，请检查求解条件\n");
            exit(1);
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
                _arb_vec_set(y_end,y_i1_ss,dim);
                break;
            }
            
            arb_add(x,x,h,prec); //改变迭代位置，这里需使用旧步长
            _arb_vec_swap(y_i1_ss,y_i,dim); //同时改变迭代初始条件
        }
        
        
        arb_mul(h,delta,h,prec); //改变下一次迭代的步长，需在X位置更新后
        
        arb_sub(t,x_end,x,prec); //用于终点判断
        
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


//在给定区间 [x_start, x_end] 中，求解 n 个点，输出这 n 个点的坐标(x,y)，每个 x 对应多个 y
//这 n 个点包括区间端点，共把区间分为 n-1 份
int ODEs_RFK45_interval_point_output(ODEs_point_output_t p_out, 
                                     my_odes_func func, const slong dim, void *param, const slong order, //常微分方程组函数
                                     const arb_t x_start, const arb_ptr y_start, //给定初始条件
                                     const arb_t x_end, //给定区间 [x_start, x_end]
                                     const slong N, const arb_t error, //N为输出点个数，区间 N-1 等分，误差为相对精度
                                     slong prec)
{
    arb_t s,D_x,x_i,x_i_p_1;
    arb_ptr v_y;
    
    arb_init(s);
    arb_init(D_x);
    arb_init(x_i);
    arb_init(x_i_p_1);
    v_y=_arb_vec_init(dim);
    
    arb_sub(D_x,x_end,x_start,prec);
    arb_div_si(D_x,D_x,N-1,prec);
    
    //这里求解时，利用了上一个点的值
    _arb_vec_set(v_y,y_start,dim);
    for(slong i=0; i<N; i++)
    {
        //x_i=x_start+(x_end-x_start)/(N-1)*i
        arb_mul_si(s,D_x,i,prec);
        arb_add(x_i,s,x_start,prec);
        
        arb_add(x_i_p_1,x_i,D_x,prec); //x_{i+1}
        
        arb_set(p_out->p_x+i,x_i);
        _arb_vec_set(p_out->p_y[i],v_y,dim);
        
        //通过[x_i,y_i],得到y_{i+1}
        ODEs_RFK45(v_y, func, dim, param, order, //常微分方程组函数
                   x_i, p_out->p_y[i],
                   x_i_p_1,2, error,
                   prec);
    }
    
    arb_clear(s);
    arb_clear(D_x);
    arb_clear(x_i);
    arb_clear(x_i_p_1);
    _arb_vec_clear(v_y,dim);
    
    return 0;
}

