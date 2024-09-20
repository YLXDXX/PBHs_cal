#include "../ODEs.h" 
#include <stdlib.h> 

//Runge–Kutta 中的 DOPRI54 方法，嵌入方法 5(4), 阶数 q=5，用于误差估算阶数 q=4  
//有 dense output，其阶数 q=4

//针对 DOPRI54 的密度输出，专门定义了一个结构体来处理

//初始化 DOPRI54 dense output 结构，为其分配内存
ODEs_DOPRI54_dense_t ODEs_DOPRI54_dense_init(slong num, slong dim)
{
    //通过 ODEs_point_output_t 声明变量后，需将该结构体指针指向结构体
    ODEs_DOPRI54_dense_t dense_out= (ODEs_DOPRI54_dense_t)calloc(1,sizeof(struct ODEs_DOPRI54_dense_output_structure));
    
    dense_out->num_keep=num;
    dense_out->dim=dim;
    
    dense_out->x=_arb_vec_init(num);
    dense_out->h=_arb_vec_init(num);
    
    dense_out->r_1=(arb_ptr*)calloc(num,sizeof(arb_ptr));
    dense_out->r_2=(arb_ptr*)calloc(num,sizeof(arb_ptr));
    dense_out->r_3=(arb_ptr*)calloc(num,sizeof(arb_ptr));
    dense_out->r_4=(arb_ptr*)calloc(num,sizeof(arb_ptr));
    dense_out->r_5=(arb_ptr*)calloc(num,sizeof(arb_ptr));
    
    for(slong i=0; i<num; i++) //对每个系数的多项式表达进行初始化
    {
        dense_out->r_1[i]=_arb_vec_init(dim);
        dense_out->r_2[i]=_arb_vec_init(dim);
        dense_out->r_3[i]=_arb_vec_init(dim);
        dense_out->r_4[i]=_arb_vec_init(dim);
        dense_out->r_5[i]=_arb_vec_init(dim);
    }
    
    return dense_out;
}

//清理 DOPRI54 dense output 结构，释放内存
void ODEs_DOPRI54_dense_clear(ODEs_DOPRI54_dense_t dense_out, slong num, slong dim)
{
    _arb_vec_clear(dense_out->x,num);
    _arb_vec_clear(dense_out->h,num);
    
    //首先释放每个多项式系数
    for(slong i=0; i<num; i++)
    {
        _arb_vec_clear(dense_out->r_1[i],dim);
        _arb_vec_clear(dense_out->r_2[i],dim);
        _arb_vec_clear(dense_out->r_3[i],dim);
        _arb_vec_clear(dense_out->r_4[i],dim);
        _arb_vec_clear(dense_out->r_5[i],dim);
    }
    
    free(dense_out);
}



// DOPRI54 计算主程序
void ODEs_DOPRI54(arb_ptr y_end, my_odes_func func, const slong dim, void *param, const slong order, //常微分方程组函数
                  const arb_t x_start, const arb_ptr y_start, //给定初始条件
                  const arb_t x_end, //求出点 x_end 对应的函数值
                  const arb_t error_abs, const arb_t error_rel, //绝对误差和相对误差
                  const ODEs_DOPRI54_dense_t dense_out, //此参数可以为空
                  slong prec)
{
    //Dormand-Prince 5(4) 方法，对应的阶数 q=5
    slong q=5;
    
    //ODEs_get_DOPRI54_cal_coe(prec); //获取系数
    //计算相关系数为全局变量，在程序运行前，在全局变量部分进行初始化
    
    arb_t gap_size; //计算区间间隔大小
    arb_init(gap_size);
    arb_sub(gap_size,x_end,x_start,prec);
    
    //计算方向
    int direction; //+1 or -1
    direction=arb_sgn_nonzero(gap_size);
    
    arb_abs(gap_size,gap_size); //需取绝对值
    
    if( arb_is_zero(gap_size) )
    {
        printf("所求点不能是初始条件点（x_end=x_start）\n");
        exit(0);
        arb_clear(gap_size);
        return;
    }
    
    
    arb_t s,t,sum,h,h_abs,x_i,x_tem,err;
    arb_init(s);
    arb_init(t);
    arb_init(sum);
    arb_init(h);
    arb_init(h_abs);
    arb_init(x_i);
    arb_init(x_tem);
    arb_init(err);
    
    arb_ptr v_s,v_t,v_w,y_i,y_i1_value,y_i1_appro;
    arb_ptr y_tem,k1,k2,k3,k4,k5,k6,k7;
    
    v_s=_arb_vec_init(dim);
    v_t=_arb_vec_init(dim);
    v_w=_arb_vec_init(dim);
    y_i=_arb_vec_init(dim);
    y_i1_value=_arb_vec_init(dim);
    y_i1_appro=_arb_vec_init(dim);
    y_tem=_arb_vec_init(dim);
    
    k1=_arb_vec_init(dim);
    k2=_arb_vec_init(dim);
    k3=_arb_vec_init(dim);
    k4=_arb_vec_init(dim);
    k5=_arb_vec_init(dim);
    k6=_arb_vec_init(dim);
    k7=_arb_vec_init(dim);
    
    
    func(v_s, x_start, y_start, dim, param, order, prec); //计算 func(x_start,y_start)，导数的初始值
    
    //计算初始迭代步长
    ODEs_select_initial_step(h_abs, func, dim, param, order, //常微分方程组函数
                             x_start, y_start, //给定初始条件
                             v_s, //func(x_start,y_start)
                             gap_size, //计算区间间隔大小 
                             q, //估算方法的阶
                             direction,
                             error_abs, error_rel, //绝对误差和相对误差
                             prec);
    
    //计算判断
    bool finish=false;
    bool dense_is_out=true;
    
    if(dense_out==NULL) //空指针不进行 dense 输出
    {
        dense_is_out=false;
    }
    
    //初始条件
    arb_set(x_i,x_start);
    _arb_vec_set(y_i,y_start,dim);
    arb_mul_si(h,h_abs,direction,prec);//h = h_abs * direction
    
    slong i=0;
    
    while( 1 )
    {
        //k1
        func(k1, x_i, y_i, dim, param, order, prec);
        _arb_vec_scalar_mul(k1,k1,dim,h,prec);
        
        //k2
        arb_mul(x_tem,ODEs_DOPRI54_coe_c+1,h,prec);
        arb_add(x_tem,x_tem,x_i,prec);
        
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOPRI54_coe_a_1+0, prec);
        _arb_vec_add(y_tem, y_i, v_s, dim, prec);
        
        func(k2, x_tem, y_tem, dim, param, order, prec);
        _arb_vec_scalar_mul(k2, k2, dim, h, prec);
        
        //k_3
        arb_mul(x_tem,ODEs_DOPRI54_coe_c+2,h,prec);
        arb_add(x_tem,x_tem,x_i,prec);
        
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOPRI54_coe_a_2+0, prec);
        _arb_vec_add(y_tem, y_i, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k2, dim, ODEs_DOPRI54_coe_a_2+1, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        func(k3, x_tem, y_tem, dim, param, order, prec);
        _arb_vec_scalar_mul(k3, k3, dim, h, prec);
        
        //k_4
        arb_mul(x_tem,ODEs_DOPRI54_coe_c+3,h,prec);
        arb_add(x_tem,x_tem,x_i,prec);
        
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOPRI54_coe_a_3+0, prec);
        _arb_vec_add(y_tem, y_i, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k2, dim, ODEs_DOPRI54_coe_a_3+1, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k3, dim, ODEs_DOPRI54_coe_a_3+2, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        func(k4, x_tem, y_tem, dim, param, order, prec);
        _arb_vec_scalar_mul(k4, k4, dim, h, prec);
        
        //k5
        arb_mul(x_tem,ODEs_DOPRI54_coe_c+4,h,prec);
        arb_add(x_tem,x_tem,x_i,prec);
        
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOPRI54_coe_a_4+0, prec);
        _arb_vec_add(y_tem, y_i, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k2, dim, ODEs_DOPRI54_coe_a_4+1, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k3, dim, ODEs_DOPRI54_coe_a_4+2, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k4, dim, ODEs_DOPRI54_coe_a_4+3, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        func(k5, x_tem, y_tem, dim, param, order, prec);
        _arb_vec_scalar_mul(k5, k5, dim, h, prec);
        
        //k6
        arb_mul(x_tem,ODEs_DOPRI54_coe_c+5,h,prec);
        arb_add(x_tem,x_tem,x_i,prec);
        
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOPRI54_coe_a_5+0, prec);
        _arb_vec_add(y_tem, y_i, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k2, dim, ODEs_DOPRI54_coe_a_5+1, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k3, dim, ODEs_DOPRI54_coe_a_5+2, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k4, dim, ODEs_DOPRI54_coe_a_5+3, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k5, dim, ODEs_DOPRI54_coe_a_5+4, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        func(k6, x_tem, y_tem, dim, param, order, prec);
        _arb_vec_scalar_mul(k6, k6, dim, h, prec);
        
        //k7
        arb_mul(x_tem,ODEs_DOPRI54_coe_c+6,h,prec);
        arb_add(x_tem,x_tem,x_i,prec);
        
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOPRI54_coe_a_6+0, prec);
        _arb_vec_add(y_tem, y_i, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k3, dim, ODEs_DOPRI54_coe_a_6+2, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k4, dim, ODEs_DOPRI54_coe_a_6+3, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k5, dim, ODEs_DOPRI54_coe_a_6+4, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k6, dim, ODEs_DOPRI54_coe_a_6+5, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        func(k7, x_tem, y_tem, dim, param, order, prec);
        _arb_vec_scalar_mul(k7, k7, dim, h, prec);
        
        
        //利用上面的 k_i 得到 y_{i+1}
        
        //y_i1_value
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOPRI54_coe_b_up+0, prec);
        _arb_vec_add(y_tem, y_i, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k3, dim, ODEs_DOPRI54_coe_b_up+2, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k4, dim, ODEs_DOPRI54_coe_b_up+3, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k5, dim, ODEs_DOPRI54_coe_b_up+4, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k6, dim, ODEs_DOPRI54_coe_b_up+5, prec);
        _arb_vec_add(y_i1_value, y_tem, v_s, dim, prec);
        
        
        
        //y_i1_appro
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOPRI54_coe_b_dw+0, prec);
        _arb_vec_add(y_tem, y_i, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k3, dim, ODEs_DOPRI54_coe_b_dw+2, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k4, dim, ODEs_DOPRI54_coe_b_dw+3, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k5, dim, ODEs_DOPRI54_coe_b_dw+4, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k6, dim, ODEs_DOPRI54_coe_b_dw+5, prec);
        _arb_vec_add(y_tem, y_tem, v_s, dim, prec);
        
        _arb_vec_scalar_mul(v_s, k7, dim, ODEs_DOPRI54_coe_b_dw+6, prec);
        _arb_vec_add(y_i1_appro, y_tem, v_s, dim, prec);
        
        
        if(finish==true)
        {
            // density output 用于插值，其精度为 h^{p-1}
            if(dense_is_out==true) //退出前保存 dense 输出相关信息
            {
                if( i+2 > dense_out->num_keep ) //这里加 2 保证后面插值时 +1 能正确执行
                {
                    printf("所分配的点个数不足，不能存入 ODEs_DOPRI54_dense_t\n请增加 num 值后重试\n");
                    exit(0);
                }
                
                //参见 https://math.stackexchange.com/questions/2947231/how-can-i-derive-the-dense-output-of-ode45
                //x_i,h_i
                arb_set(dense_out->x+i,x_i);
                arb_set(dense_out->h+i,h);
                
                //计算 r_1 到 r_5 存入结构体
                //r_1 = y_n
                _arb_vec_set(dense_out->r_1[i],y_i,dim);
                
                //r2=yn+1−yn
                _arb_vec_sub(dense_out->r_2[i],y_i1_value,y_i,dim,prec);
                
                
                //r3=yn+hfn−yn+1，其中 hfn 为 k1
                _arb_vec_add(v_s,y_i,k1,dim,prec);
                _arb_vec_sub(dense_out->r_3[i],v_s,y_i1_value,dim,prec);
                
                
                //r4=2(yn+1−yn)−h(fn+fn+1) 其中 f_{n+1}(x_{i+1}, y_{i+1})
                arb_add(s,x_i,h,prec);
                func(v_w, s, y_i1_value, dim, param, order, prec);
                _arb_vec_scalar_mul(v_w,v_w,dim,h,prec);
                _arb_vec_add(v_s,v_w,k1,dim,prec);
                
                arb_set_si(t,2);
                _arb_vec_scalar_mul(v_t,dense_out->r_2[i],dim,t,prec);
                _arb_vec_sub(dense_out->r_4[i],v_t,v_s,dim,prec);
                
                
                //r5=d1hfn+d3k3+d4k4+d5k5+d6k6+d7hfn+1
                _arb_vec_scalar_mul(v_s,k1,dim,ODEs_DOPRI54_coe_dense_out+0,prec);
                _arb_vec_scalar_mul(v_t,k3,dim,ODEs_DOPRI54_coe_dense_out+2,prec);
                _arb_vec_add(v_s,v_s,v_t,dim,prec);
                
                _arb_vec_scalar_mul(v_t,k4,dim,ODEs_DOPRI54_coe_dense_out+3,prec);
                _arb_vec_add(v_s,v_s,v_t,dim,prec);
                
                _arb_vec_scalar_mul(v_t,k5,dim,ODEs_DOPRI54_coe_dense_out+4,prec);
                _arb_vec_add(v_s,v_s,v_t,dim,prec);
                
                _arb_vec_scalar_mul(v_t,k6,dim,ODEs_DOPRI54_coe_dense_out+5,prec);
                _arb_vec_add(v_s,v_s,v_t,dim,prec);
                
                _arb_vec_scalar_mul(v_t,v_w,dim,ODEs_DOPRI54_coe_dense_out+6,prec);
                _arb_vec_add(dense_out->r_5[i],v_s,v_t,dim,prec);
                
                //需要对最后一个点 x_end 进行特殊处理，方便后面查找 x_end=x+2*h
                //最后一点只需存 x 即可
                arb_mul_si(s,h,2,prec);
                arb_add(dense_out->x+i+1,x_i,s,prec);
                
                dense_out->num_real=i+2; //记录真实存入点的个数
                
                //i++; //更新索引
            }
            
            //完成，退出
            //arb_set(x_i,x_end);
            _arb_vec_set(y_end,y_i1_value,dim);
            break;
        }
        
        //再利用上面的值，判断误差
        _arb_vec_sub(v_s,y_i1_value,y_i1_appro,dim,prec);
        
        arb_zero(sum); 
        for(slong i=0; i<dim; i++)
        {
            //sci = Atoli + max(|y0i|, |y1i|) · Rtoli
            arb_abs(s,y_i+i);
            arb_abs(t,y_i1_value+i);
            arb_max(s,s,t,prec);
            arb_mul(s,s,error_rel,prec);
            arb_add(s,s,error_abs,prec);
            
            arb_div(s,v_s+i,s,prec);
            arb_sqr(s,s,prec);
            arb_add(sum,sum,s,prec);
            
            arb_get_mid_arb(y_i1_value+i,y_i1_value+i); //迭代次数足够多时，计算误差会很大，消去
        }
        arb_div_si(s,sum,dim,prec);
        arb_sqrt(err,s,prec);
        
        
        if(!arb_is_finite(err))
        {
            printf("遇到错误已退出，可能原因：\n");
            printf("1) 精度要求太高，当前精度不能满足要求，请提高 prec 值后重试\n");
            printf("2) 求解条件设定有问题，请检查求解条件\n");
            exit(1);
        }
        
        arb_one(s);
        if( arb_le(err,s) ) // err≤1 满足误差要求，进入下一步
        {
            // density output 用于插值，其精度为 h^{p-1}
            if(dense_is_out==true)
            {
                if( i+2 > dense_out->num_keep ) //这里加 2 保证后面插值时 +1 能正确执行
                {
                    printf("所分配的点个数不足，不能存入 ODEs_DOPRI54_dense_t\n请增加 num 值后重试\n");
                    exit(0);
                }
                
                //参见 https://math.stackexchange.com/questions/2947231/how-can-i-derive-the-dense-output-of-ode45
                //x_i,h_i
                arb_set(dense_out->x+i,x_i);
                arb_set(dense_out->h+i,h);
                
                //计算 r_1 到 r_5 存入结构体
                //r_1 = y_n
                _arb_vec_set(dense_out->r_1[i],y_i,dim);
                
                //r2=yn+1−yn
                _arb_vec_sub(dense_out->r_2[i],y_i1_value,y_i,dim,prec);
                
                
                //r3=yn+hfn−yn+1，其中 hfn 为 k1
                _arb_vec_add(v_s,y_i,k1,dim,prec);
                _arb_vec_sub(dense_out->r_3[i],v_s,y_i1_value,dim,prec);
                
                
                //r4=2(yn+1−yn)−h(fn+fn+1) 其中 f_{n+1}(x_{i+1}, y_{i+1})
                arb_add(s,x_i,h,prec);
                func(v_w, s, y_i1_value, dim, param, order, prec);
                _arb_vec_scalar_mul(v_w,v_w,dim,h,prec);
                _arb_vec_add(v_s,v_w,k1,dim,prec);
                
                arb_set_si(t,2);
                _arb_vec_scalar_mul(v_t,dense_out->r_2[i],dim,t,prec);
                _arb_vec_sub(dense_out->r_4[i],v_t,v_s,dim,prec);
                
                
                //r5=d1hfn+d3k3+d4k4+d5k5+d6k6+d7hfn+1
                _arb_vec_scalar_mul(v_s,k1,dim,ODEs_DOPRI54_coe_dense_out+0,prec);
                _arb_vec_scalar_mul(v_t,k3,dim,ODEs_DOPRI54_coe_dense_out+2,prec);
                _arb_vec_add(v_s,v_s,v_t,dim,prec);
                
                _arb_vec_scalar_mul(v_t,k4,dim,ODEs_DOPRI54_coe_dense_out+3,prec);
                _arb_vec_add(v_s,v_s,v_t,dim,prec);
                
                _arb_vec_scalar_mul(v_t,k5,dim,ODEs_DOPRI54_coe_dense_out+4,prec);
                _arb_vec_add(v_s,v_s,v_t,dim,prec);
                
                _arb_vec_scalar_mul(v_t,k6,dim,ODEs_DOPRI54_coe_dense_out+5,prec);
                _arb_vec_add(v_s,v_s,v_t,dim,prec);
                
                _arb_vec_scalar_mul(v_t,v_w,dim,ODEs_DOPRI54_coe_dense_out+6,prec);
                _arb_vec_add(dense_out->r_5[i],v_s,v_t,dim,prec);
                
                i++; //更新索引
            }
            
            ///更新迭代位置
            arb_add(x_i,x_i,h,prec); //这里的 h 需要用老位置
            _arb_vec_swap(y_i1_value,y_i,dim); //同时改变迭代初始条件
            
            //预估下一迭代的步长
            arb_inv(s,err,prec); //求新步长，见 (4.12)
            arb_root_ui(s,s,q+1,prec);
            arb_mul(s,s,ODEs_step_factor,prec);
            
            arb_max(s,s,ODEs_step_min_factor,prec);
            arb_min(s,s,ODEs_step_max_factor,prec);
            
            arb_mul(h_abs,h_abs,s,prec); //新步长
            arb_get_mid_arb(h_abs,h_abs);
            
            //利用预估的步长和 |x_end-x_i| 判断，是否终止
            arb_sub(s,x_end,x_i,prec);
            arb_abs(s,s);
            
            if( arb_le(s,h_abs) )
            {
                //在下一次迭代后，将退出，完成计算
                arb_set(h_abs,s);
                finish=true;
            }
            
            arb_mul_si(h,h_abs,direction,prec);//h = h_abs * direction
            
        }else
        {
            //没达到精度，在 x_i 停留，减小步长，继续
            //同样，预估下一迭代的步长
            arb_inv(s,err,prec); //求新步长，见 (4.12)
            arb_root_ui(s,s,q+1,prec);
            arb_mul(s,s,ODEs_step_factor,prec);
            
            arb_max(s,s,ODEs_step_min_factor,prec);
            arb_min(s,s,ODEs_step_max_factor,prec);
            
            arb_mul(h_abs,h_abs,s,prec); //新步长
            arb_get_mid_arb(h_abs,h_abs);
            
            arb_mul_si(h,h_abs,direction,prec);//h = h_abs * direction
        }
    }
    
    arb_clear(gap_size);
    arb_clear(s);
    arb_clear(t);
    arb_clear(sum);
    arb_clear(h);
    arb_clear(h_abs);
    arb_clear(x_i);
    arb_clear(x_tem);
    arb_clear(err);
    
    _arb_vec_clear(v_s,dim);
    _arb_vec_clear(v_t,dim);
    _arb_vec_clear(v_w,dim);
    _arb_vec_clear(y_i,dim);
    _arb_vec_clear(y_i1_value,dim);
    _arb_vec_clear(y_i1_appro,dim);
    _arb_vec_clear(y_tem,dim);
    
    _arb_vec_clear(k1,dim);
    _arb_vec_clear(k2,dim);
    _arb_vec_clear(k3,dim);
    _arb_vec_clear(k4,dim);
    _arb_vec_clear(k5,dim);
    _arb_vec_clear(k6,dim);
    _arb_vec_clear(k7,dim);
    
}

