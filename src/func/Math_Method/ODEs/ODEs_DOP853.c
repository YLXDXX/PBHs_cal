#include "../ODEs.h" 
#include <stdlib.h> 

//Runge–Kutta 中的 DOP853 方法，嵌入方法 8(5,3), 阶数 q=8，用于误差估算阶数 q=3/5  
//有 dense output，其阶数 q=7

//注意，此方法受系数求解精度的限制
//当计算的有效数位达到 30 位时，精度将不能再提高（与系数的精度保持一致）


//针对 DOP853 的密度输出，专门定义了一个结构体来处理

//初始化 DOP853 dense output 结构，为其分配内存
ODEs_DOP853_dense_t ODEs_DOP853_dense_init(slong num, slong dim)
{
    //通过 ODEs_point_output_t 声明变量后，需将该结构体指针指向结构体
    ODEs_DOP853_dense_t dense_out= (ODEs_DOP853_dense_t)calloc(1,sizeof(struct ODEs_DOP853_dense_output_structure));
    
    dense_out->num_keep=num;
    dense_out->dim=dim;
    
    dense_out->x=_arb_vec_init(num);
    dense_out->h=_arb_vec_init(num);
    
    dense_out->rcont1=(arb_ptr*)calloc(num,sizeof(arb_ptr));
    dense_out->rcont2=(arb_ptr*)calloc(num,sizeof(arb_ptr));
    dense_out->rcont3=(arb_ptr*)calloc(num,sizeof(arb_ptr));
    dense_out->rcont4=(arb_ptr*)calloc(num,sizeof(arb_ptr));
    dense_out->rcont5=(arb_ptr*)calloc(num,sizeof(arb_ptr));
    dense_out->rcont6=(arb_ptr*)calloc(num,sizeof(arb_ptr));
    dense_out->rcont7=(arb_ptr*)calloc(num,sizeof(arb_ptr));
    dense_out->rcont8=(arb_ptr*)calloc(num,sizeof(arb_ptr));
    
    for(slong i=0; i<num; i++) //对每个系数的多项式表达进行初始化
    {
        dense_out->rcont1[i]=_arb_vec_init(dim);
        dense_out->rcont2[i]=_arb_vec_init(dim);
        dense_out->rcont3[i]=_arb_vec_init(dim);
        dense_out->rcont4[i]=_arb_vec_init(dim);
        dense_out->rcont5[i]=_arb_vec_init(dim);
        dense_out->rcont6[i]=_arb_vec_init(dim);
        dense_out->rcont7[i]=_arb_vec_init(dim);
        dense_out->rcont8[i]=_arb_vec_init(dim);
    }
    
    return dense_out;
}

//清理 DOP853 dense output 结构，释放内存
void ODEs_DOP853_dense_clear(ODEs_DOP853_dense_t dense_out, slong num, slong dim)
{
    _arb_vec_clear(dense_out->x,num);
    _arb_vec_clear(dense_out->h,num);
    
    //首先释放每个多项式系数
    for(slong i=0; i<num; i++)
    {
        _arb_vec_clear(dense_out->rcont1[i],dim);
        _arb_vec_clear(dense_out->rcont2[i],dim);
        _arb_vec_clear(dense_out->rcont3[i],dim);
        _arb_vec_clear(dense_out->rcont4[i],dim);
        _arb_vec_clear(dense_out->rcont5[i],dim);
        _arb_vec_clear(dense_out->rcont6[i],dim);
        _arb_vec_clear(dense_out->rcont7[i],dim);
        _arb_vec_clear(dense_out->rcont8[i],dim);
    }
    
    free(dense_out);
}



// DOP853 计算主程序
//此程序参见 https://github.com/robclewley/pydstool/blob/master/PyDSTool/integrator/dop853.c
void ODEs_DOP853(arb_ptr y_end, my_odes_func func, const slong dim, void *param, const slong order, //常微分方程组函数
                  const arb_t x_start, const arb_ptr y_start, //给定初始条件
                  const arb_t x_end, //求出点 x_end 对应的函数值
                  const arb_t error_abs, const arb_t error_rel, //绝对误差和相对误差
                  const ODEs_DOP853_dense_t dense_out, //此参数可以为空
                  slong prec)
{
    //Hairer's 8/5/3 adaption of the Dormand-Prince Runge-Kutta method. (7th order interpolant).
    slong q=8;
    
    //ODEs_get_DOP853_cal_coe(prec); //获取系数
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
    
    
    arb_t s,t,sci,sum,sum_2,h,h_abs,x_i,x_tem,err;
    arb_init(s);
    arb_init(t);
    arb_init(sci);
    arb_init(sum);
    arb_init(sum_2);
    arb_init(h);
    arb_init(h_abs);
    arb_init(x_i);
    arb_init(x_tem);
    arb_init(err);
    
    arb_ptr v_s,v_t,v_w,y_i,y_i1_value,y_i1_diff_3,y_i1_diff_5;
    arb_ptr y_tem,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k14,k15,k16;
    
    v_s=_arb_vec_init(dim);
    v_t=_arb_vec_init(dim);
    v_w=_arb_vec_init(dim);
    y_i=_arb_vec_init(dim);
    y_i1_value=_arb_vec_init(dim);
    y_i1_diff_3=_arb_vec_init(dim);
    y_i1_diff_5=_arb_vec_init(dim);
    y_tem=_arb_vec_init(dim);
    
    k1=_arb_vec_init(dim);
    k2=_arb_vec_init(dim);
    k3=_arb_vec_init(dim);
    k4=_arb_vec_init(dim);
    k5=_arb_vec_init(dim);
    k6=_arb_vec_init(dim);
    k7=_arb_vec_init(dim);
    k8=_arb_vec_init(dim);
    k9=_arb_vec_init(dim);
    k10=_arb_vec_init(dim);
    k11=_arb_vec_init(dim);
    k12=_arb_vec_init(dim);
    k14=_arb_vec_init(dim);
    k15=_arb_vec_init(dim);
    k16=_arb_vec_init(dim);
    
    
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
        //这里的 k_i 与其它方法不同，并没有乘以 h
        //k1
        func(k1, x_i, y_i, dim, param, order, prec);
        
        //k2
        arb_mul(x_tem,ODEs_DOP853_coe_c2,h,prec);
        arb_add(x_tem,x_tem,x_i,prec);
        
        //y[i] + h * a21 * k1[i]
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_a0201, prec);
        
        _arb_vec_scalar_mul(v_s, v_s, dim, h, prec);
        _arb_vec_add(y_tem, y_i, v_s, dim, prec);
        
        func(k2, x_tem, y_tem, dim, param, order, prec);
        
        //k_3
        arb_mul(x_tem,ODEs_DOP853_coe_c3,h,prec);
        arb_add(x_tem,x_tem,x_i,prec);
        
        //y[i] + h * (a31*k1[i] + a32*k2[i])
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_a0301, prec);
        _arb_vec_scalar_mul(v_w, k2, dim, ODEs_DOP853_coe_a0302, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        
        _arb_vec_scalar_mul(v_s, v_s, dim, h, prec);
        _arb_vec_add(y_tem, y_i, v_s, dim, prec);
        
        func(k3, x_tem, y_tem, dim, param, order, prec);
        
        //k_4
        arb_mul(x_tem,ODEs_DOP853_coe_c4,h,prec);
        arb_add(x_tem,x_tem,x_i,prec);
        
        //y[i] + h * (a41*k1[i] + a43*k3[i])
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_a0401, prec);
        _arb_vec_scalar_mul(v_w, k3, dim, ODEs_DOP853_coe_a0403, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        
        _arb_vec_scalar_mul(v_s, v_s, dim, h, prec);
        _arb_vec_add(y_tem, y_i, v_s, dim, prec);
        
        func(k4, x_tem, y_tem, dim, param, order, prec);
        
        //k5
        arb_mul(x_tem,ODEs_DOP853_coe_c5,h,prec);
        arb_add(x_tem,x_tem,x_i,prec);
        
        //y[i] + h * (a51*k1[i] + a53*k3[i] + a54*k4[i])
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_a0501, prec);
        _arb_vec_scalar_mul(v_w, k3, dim, ODEs_DOP853_coe_a0503, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k4, dim, ODEs_DOP853_coe_a0504, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        
        _arb_vec_scalar_mul(v_s, v_s, dim, h, prec);
        _arb_vec_add(y_tem, y_i, v_s, dim, prec);
        
        func(k5, x_tem, y_tem, dim, param, order, prec);
        
        //k6
        arb_mul(x_tem,ODEs_DOP853_coe_c6,h,prec);
        arb_add(x_tem,x_tem,x_i,prec);
        
        //y[i] + h * (a61*k1[i] + a64*k4[i] + a65*k5[i])
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_a0601, prec);
        _arb_vec_scalar_mul(v_w, k4, dim, ODEs_DOP853_coe_a0604, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k5, dim, ODEs_DOP853_coe_a0605, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        
        _arb_vec_scalar_mul(v_s, v_s, dim, h, prec);
        _arb_vec_add(y_tem, y_i, v_s, dim, prec);
        
        func(k6, x_tem, y_tem, dim, param, order, prec);
        
        //k7
        arb_mul(x_tem,ODEs_DOP853_coe_c7,h,prec);
        arb_add(x_tem,x_tem,x_i,prec);
        
        //y[i] + h * (a71*k1[i] + a74*k4[i] + a75*k5[i] + a76*k6[i])
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_a0701, prec);
        _arb_vec_scalar_mul(v_w, k4, dim, ODEs_DOP853_coe_a0704, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k5, dim, ODEs_DOP853_coe_a0705, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_a0706, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        
        _arb_vec_scalar_mul(v_s, v_s, dim, h, prec);
        _arb_vec_add(y_tem, y_i, v_s, dim, prec);
        
        func(k7, x_tem, y_tem, dim, param, order, prec);
        
        //k8
        arb_mul(x_tem,ODEs_DOP853_coe_c8,h,prec);
        arb_add(x_tem,x_tem,x_i,prec);
        
        //y[i] + h * (a81*k1[i] + a84*k4[i] + a85*k5[i] + a86*k6[i] + a87*k7[i])
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_a0801, prec);
        _arb_vec_scalar_mul(v_w, k4, dim, ODEs_DOP853_coe_a0804, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k5, dim, ODEs_DOP853_coe_a0805, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_a0806, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_a0807, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        
        _arb_vec_scalar_mul(v_s, v_s, dim, h, prec);
        _arb_vec_add(y_tem, y_i, v_s, dim, prec);
        
        func(k8, x_tem, y_tem, dim, param, order, prec);
        
        //k9
        arb_mul(x_tem,ODEs_DOP853_coe_c9,h,prec);
        arb_add(x_tem,x_tem,x_i,prec);
        
        //y[i] + h * (a91*k1[i] + a94*k4[i] + a95*k5[i] + a96*k6[i] + a97*k7[i] + a98*k8[i])
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_a0901, prec);
        _arb_vec_scalar_mul(v_w, k4, dim, ODEs_DOP853_coe_a0904, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k5, dim, ODEs_DOP853_coe_a0905, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_a0906, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_a0907, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_a0908, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        
        _arb_vec_scalar_mul(v_s, v_s, dim, h, prec);
        _arb_vec_add(y_tem, y_i, v_s, dim, prec);
        
        func(k9, x_tem, y_tem, dim, param, order, prec);
        
        //k10
        arb_mul(x_tem,ODEs_DOP853_coe_c10,h,prec);
        arb_add(x_tem,x_tem,x_i,prec);
        
        //y[i] + h * (a101*k1[i] + a104*k4[i] + a105*k5[i] + a106*k6[i] + a107*k7[i] + a108*k8[i] + a109*k9[i])
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_a1001, prec);
        _arb_vec_scalar_mul(v_w, k4, dim, ODEs_DOP853_coe_a1004, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k5, dim, ODEs_DOP853_coe_a1005, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_a1006, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_a1007, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_a1008, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k9, dim, ODEs_DOP853_coe_a1009, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        
        _arb_vec_scalar_mul(v_s, v_s, dim, h, prec);
        _arb_vec_add(y_tem, y_i, v_s, dim, prec);
        
        func(k10, x_tem, y_tem, dim, param, order, prec);
        
        //原 c 程序为了节约内存，重复利用了后面用不到的 k2 k3 k4 k5
        //这里为了清晰，不再采用
        //k11 (k2)
        arb_mul(x_tem,ODEs_DOP853_coe_c11,h,prec);
        arb_add(x_tem,x_tem,x_i,prec);
        
        //y[i] + h * (a111*k1[i] + a114*k4[i] + a115*k5[i] + a116*k6[i] + a117*k7[i] + a118*k8[i] + a119*k9[i] + a1110*k10[i])
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_a1101, prec);
        _arb_vec_scalar_mul(v_w, k4, dim, ODEs_DOP853_coe_a1104, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k5, dim, ODEs_DOP853_coe_a1105, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_a1106, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_a1107, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_a1108, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k9, dim, ODEs_DOP853_coe_a1109, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k10, dim, ODEs_DOP853_coe_a1110, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        
        _arb_vec_scalar_mul(v_s, v_s, dim, h, prec);
        _arb_vec_add(y_tem, y_i, v_s, dim, prec);
        
        func(k11, x_tem, y_tem, dim, param, order, prec);
        
        //k12 (k3)
        arb_add(x_tem,x_i,h,prec); // x_i +h
        
        //y[i] + h * (a121*k1[i] + a124*k4[i] + a125*k5[i] + a126*k6[i] + a127*k7[i] + a128*k8[i] + a129*k9[i] +
        //            a1210*k10[i] + a1211*k2[i])
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_a1201, prec);
        _arb_vec_scalar_mul(v_w, k4, dim, ODEs_DOP853_coe_a1204, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k5, dim, ODEs_DOP853_coe_a1205, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_a1206, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_a1207, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_a1208, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k9, dim, ODEs_DOP853_coe_a1209, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k10, dim, ODEs_DOP853_coe_a1210, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k11, dim, ODEs_DOP853_coe_a1211, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        
        _arb_vec_scalar_mul(v_s, v_s, dim, h, prec);
        _arb_vec_add(y_tem, y_i, v_s, dim, prec);
        
        func(k12, x_tem, y_tem, dim, param, order, prec);
        
        //将 k4 利用起来，后面复用
        //k4
        //k4[i] = b1*k1[i] + b6*k6[i] + b7*k7[i] + b8*k8[i] + b9*k9[i] +b10*k10[i] + b11*k2[i] + b12*k3[i];
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_b1, prec);
        _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_b6, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_b7, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_b8, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k9, dim, ODEs_DOP853_coe_b9, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k10, dim, ODEs_DOP853_coe_b10, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k11, dim, ODEs_DOP853_coe_b11, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k12, dim, ODEs_DOP853_coe_b12, prec);
        _arb_vec_add(k4,v_s,v_w,dim,prec);
        
        
        //利用前面的值，得到近似值 y_i1_value (k5)
        //k5[i] = y[i] + h * k4[i], 即最后得到的 y_{i+1} 值
        _arb_vec_scalar_mul(v_s, k4, dim, h, prec);
        _arb_vec_add(y_i1_value, y_i, v_s, dim, prec);
        
        
        if(finish==true)
        {
            // density output 用于插值，其精度为 h^{p-1}
            if(dense_is_out==true) //退出前保存 dense 输出相关信息
            {
                if( i+2 > dense_out->num_keep ) //这里加 2 保证后面插值时 +1 能正确执行
                {
                    printf("所分配的点个数不足，不能存入 ODEs_DOP853_dense_t\n请增加 num 值后重试\n");
                    exit(0);
                }
                
                //x_i,h_i
                arb_set(dense_out->x+i,x_i);
                arb_set(dense_out->h+i,h);
                
                //计算 rcont1 到 rcont8 存入结构体
                //rcont1[i] = y[i];
                _arb_vec_set(dense_out->rcont1[i],y_i,dim);
                
                //rcont2[i] = ydiff = k5[i] - y[i]
                _arb_vec_sub(dense_out->rcont2[i],y_i1_value,y_i,dim,prec);
                
                //rcont3[i] = bspl = h * k1[i] - ydiff;
                _arb_vec_scalar_mul(v_s,k1,dim,h,prec);
                _arb_vec_sub(dense_out->rcont3[i],v_s,dense_out->rcont2[i],dim,prec);
                
                //rcont4[i] = ydiff - h*k4[i] - bspl;
                arb_add(s,x_i,h,prec); //fcn (n, xph, k5, pars, k4);
                func(k4, s, y_i1_value, dim, param, order, prec); //对 k4 进行更新，相当于下一次迭代的 k1
                
                _arb_vec_scalar_mul(v_s,k4,dim,h,prec);
                _arb_vec_sub(v_s,dense_out->rcont2[i],v_s,dim,prec);
                _arb_vec_sub(dense_out->rcont4[i],v_s,dense_out->rcont3[i],dim,prec);
                
                //rcont5[i] = d41*k1[i] + d46*k6[i] + d47*k7[i] + d48*k8[i] + 
                //            d49*k9[i] + d410*k10[i] + d411*k2[i] + d412*k3[i];
                _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_d401, prec);
                _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_d406, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_d407, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_d408, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k9, dim, ODEs_DOP853_coe_d409, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k10, dim, ODEs_DOP853_coe_d410, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k11, dim, ODEs_DOP853_coe_d411, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k12, dim, ODEs_DOP853_coe_d412, prec);
                _arb_vec_add(dense_out->rcont5[i],v_s,v_w,dim,prec);
                
                //rcont6[i] = d51*k1[i] + d56*k6[i] + d57*k7[i] + d58*k8[i] +
                //            d59*k9[i] + d510*k10[i] + d511*k2[i] + d512*k3[i];
                _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_d501, prec);
                _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_d506, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_d507, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_d508, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k9, dim, ODEs_DOP853_coe_d509, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k10, dim, ODEs_DOP853_coe_d510, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k11, dim, ODEs_DOP853_coe_d511, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k12, dim, ODEs_DOP853_coe_d512, prec);
                _arb_vec_add(dense_out->rcont6[i],v_s,v_w,dim,prec);
                
                
                //rcont7[i] = d61*k1[i] + d66*k6[i] + d67*k7[i] + d68*k8[i] +
                //            d69*k9[i] + d610*k10[i] + d611*k2[i] + d612*k3[i];
                _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_d601, prec);
                _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_d606, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_d607, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_d608, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k9, dim, ODEs_DOP853_coe_d609, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k10, dim, ODEs_DOP853_coe_d610, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k11, dim, ODEs_DOP853_coe_d611, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k12, dim, ODEs_DOP853_coe_d612, prec);
                _arb_vec_add(dense_out->rcont7[i],v_s,v_w,dim,prec);
                
                
                //rcont8[i] = d71*k1[i] + d76*k6[i] + d77*k7[i] + d78*k8[i] +
                //            d79*k9[i] + d710*k10[i] + d711*k2[i] + d712*k3[i];
                _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_d701, prec);
                _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_d706, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_d707, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_d708, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k9, dim, ODEs_DOP853_coe_d709, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k10, dim, ODEs_DOP853_coe_d710, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k11, dim, ODEs_DOP853_coe_d711, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k12, dim, ODEs_DOP853_coe_d712, prec);
                _arb_vec_add(dense_out->rcont7[i],v_s,v_w,dim,prec);
                
                
                /* the next three function evaluations */
                //k14 (k10)
                arb_mul(x_tem,ODEs_DOP853_coe_c14,h,prec);
                arb_add(x_tem,x_tem,x_i,prec);
                
                //y[i] + h * (a141*k1[i] + a147*k7[i] + a148*k8[i] + a149*k9[i] + a1410*k10[i] + a1411*k2[i] + 
                //            a1412*k3[i] + a1413*k4[i])
                _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_a1401, prec);
                _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_a1407, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_a1408, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k9, dim, ODEs_DOP853_coe_a1409, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k10, dim, ODEs_DOP853_coe_a1410, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k11, dim, ODEs_DOP853_coe_a1411, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k12, dim, ODEs_DOP853_coe_a1412, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k4, dim, ODEs_DOP853_coe_a1413, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                
                _arb_vec_scalar_mul(v_s, v_s, dim, h, prec);
                _arb_vec_add(y_tem, y_i, v_s, dim, prec);
                
                func(k14, x_tem, y_tem, dim, param, order, prec);
                
                //k15 (k2)
                arb_mul(x_tem,ODEs_DOP853_coe_c15,h,prec);
                arb_add(x_tem,x_tem,x_i,prec);
                
                //y[i] + h * (a151*k1[i] + a156*k6[i] + a157*k7[i] + a158*k8[i] +
                //            a1511*k2[i] + a1512*k3[i] + a1513*k4[i] + a1514*k10[i]);
                _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_a1501, prec);
                _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_a1506, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_a1507, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_a1508, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k11, dim, ODEs_DOP853_coe_a1511, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k12, dim, ODEs_DOP853_coe_a1512, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k4, dim, ODEs_DOP853_coe_a1513, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k14, dim, ODEs_DOP853_coe_a1514, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                
                _arb_vec_scalar_mul(v_s, v_s, dim, h, prec);
                _arb_vec_add(y_tem, y_i, v_s, dim, prec);
                
                func(k15, x_tem, y_tem, dim, param, order, prec);
                
                //k16 (k3)
                arb_mul(x_tem,ODEs_DOP853_coe_c16,h,prec);
                arb_add(x_tem,x_tem,x_i,prec);
                
                //y[i] + h * (a161*k1[i] + a166*k6[i] + a167*k7[i] + a168*k8[i] +
                //            a169*k9[i] + a1613*k4[i] + a1614*k10[i] + a1615*k2[i]);
                _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_a1601, prec);
                _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_a1606, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_a1607, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_a1608, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k9, dim, ODEs_DOP853_coe_a1609, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k4, dim, ODEs_DOP853_coe_a1613, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k14, dim, ODEs_DOP853_coe_a1614, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k15, dim, ODEs_DOP853_coe_a1615, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                
                _arb_vec_scalar_mul(v_s, v_s, dim, h, prec);
                _arb_vec_add(y_tem, y_i, v_s, dim, prec);
                
                func(k16, x_tem, y_tem, dim, param, order, prec);
                
                //利用 k14 k15 k16 再求 dense output 系数
                
                //rcont5[i] = h * (rcont5[i] + d413*k4[i] + d414*k10[i] + d415*k2[i] + d416*k3[i]);
                _arb_vec_scalar_mul(v_s, k4, dim, ODEs_DOP853_coe_d413, prec);
                _arb_vec_scalar_mul(v_w, k14, dim, ODEs_DOP853_coe_d414, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k15, dim, ODEs_DOP853_coe_d415, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k16, dim, ODEs_DOP853_coe_d416, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                
                _arb_vec_add(v_s,v_s,dense_out->rcont5[i],dim,prec);
                _arb_vec_scalar_mul(dense_out->rcont5[i], v_s, dim, h, prec);
                
                //rcont6[i] = h * (rcont6[i] + d513*k4[i] + d514*k10[i] + d515*k2[i] + d516*k3[i]);
                _arb_vec_scalar_mul(v_s, k4, dim, ODEs_DOP853_coe_d513, prec);
                _arb_vec_scalar_mul(v_w, k14, dim, ODEs_DOP853_coe_d514, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k15, dim, ODEs_DOP853_coe_d515, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k16, dim, ODEs_DOP853_coe_d516, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                
                _arb_vec_add(v_s,v_s,dense_out->rcont6[i],dim,prec);
                _arb_vec_scalar_mul(dense_out->rcont6[i], v_s, dim, h, prec);
                
                //rcont7[i] = h * (rcont7[i] + d613*k4[i] + d614*k10[i] + d615*k2[i] + d616*k3[i]);
                _arb_vec_scalar_mul(v_s, k4, dim, ODEs_DOP853_coe_d613, prec);
                _arb_vec_scalar_mul(v_w, k14, dim, ODEs_DOP853_coe_d614, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k15, dim, ODEs_DOP853_coe_d615, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k16, dim, ODEs_DOP853_coe_d616, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                
                _arb_vec_add(v_s,v_s,dense_out->rcont7[i],dim,prec);
                _arb_vec_scalar_mul(dense_out->rcont7[i], v_s, dim, h, prec);
                
                //rcont8[i] = h * (rcont8[i] + d713*k4[i] + d714*k10[i] + d715*k2[i] + d716*k3[i]);
                _arb_vec_scalar_mul(v_s, k4, dim, ODEs_DOP853_coe_d713, prec);
                _arb_vec_scalar_mul(v_w, k14, dim, ODEs_DOP853_coe_d714, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k15, dim, ODEs_DOP853_coe_d715, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k16, dim, ODEs_DOP853_coe_d716, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                
                _arb_vec_add(v_s,v_s,dense_out->rcont8[i],dim,prec);
                _arb_vec_scalar_mul(dense_out->rcont8[i], v_s, dim, h, prec);
                
                
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
        
        //前面是利用 (y_{i+1}-y'_{i+1}) 估计误差
        //而这个差值，在这里，可直接利用系数得到
        
        //这里估算误差用了两个值，一个是三次的，一个是五次的
        //三次的差值 k4[i] - bhh1*k1[i] - bhh2*k9[i] - bhh3*k3[i]
        _arb_vec_scalar_mul(v_w,k1,dim,ODEs_DOP853_coe_bhh1,prec);
        _arb_vec_sub(v_s,k4,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w,k9,dim,ODEs_DOP853_coe_bhh2,prec);
        _arb_vec_sub(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w,k12,dim,ODEs_DOP853_coe_bhh3,prec);
        _arb_vec_sub(y_i1_diff_3,v_s,v_w,dim,prec);
        
        //五次的 er1*k1[i] + er6*k6[i] + er7*k7[i] + er8*k8[i] + er9*k9[i] + er10 * k10[i] + er11*k2[i] + er12*k3[i]
        _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_er1, prec);
        _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_er6, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_er7, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_er8, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k9, dim, ODEs_DOP853_coe_er9, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k10, dim, ODEs_DOP853_coe_er10, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k11, dim, ODEs_DOP853_coe_er11, prec);
        _arb_vec_add(v_s,v_s,v_w,dim,prec);
        _arb_vec_scalar_mul(v_w, k12, dim, ODEs_DOP853_coe_er12, prec);
        _arb_vec_add(y_i1_diff_5,v_s,v_w,dim,prec);
        
        
        arb_zero(sum);
        arb_zero(sum_2);
        for(slong i=0; i<dim; i++)
        {
            //sci = Atoli + max(|y0i|, |y1i|) · Rtoli
            arb_abs(s,y_i+i);
            arb_abs(t,y_i1_value+i);
            arb_max(s,s,t,prec);
            arb_mul(s,s,error_rel,prec);
            arb_add(sci,s,error_abs,prec);
            
            //三次的 
            arb_div(s,y_i1_diff_3+i,sci,prec);
            arb_sqr(s,s,prec);
            arb_add(sum_2,sum_2,s,prec); //err2 += sqr*sqr;
            
            //五次的
            arb_div(s,y_i1_diff_5+i,sci,prec);
            arb_sqr(s,s,prec);
            arb_add(sum,sum,s,prec); //err += sqr*sqr;
            
            arb_get_mid_arb(y_i1_value+i,y_i1_value+i); //迭代次数足够多时，计算误差会很大，消去
        }
        
        //deno = err + 0.01 * err2
        //err = fabs(h) * err * sqrt( 1.0 / (deno*dim) );
        arb_div_ui(s,sum_2,100,prec);
        arb_add(s,s,sum,prec);
        arb_mul_si(s,s,dim,prec);
        arb_sqrt(s,s,prec);
        
        arb_div(s,sum,s,prec);
        arb_mul(err,s,h_abs,prec);
        
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
                    printf("所分配的点个数不足，不能存入 ODEs_DOP853_dense_t\n请增加 num 值后重试\n");
                    exit(0);
                }
                
                //x_i,h_i
                arb_set(dense_out->x+i,x_i);
                arb_set(dense_out->h+i,h);
                
                //计算 rcont1 到 rcont8 存入结构体
                //rcont1[i] = y[i];
                _arb_vec_set(dense_out->rcont1[i],y_i,dim);
                
                //rcont2[i] = ydiff = k5[i] - y[i]
                _arb_vec_sub(dense_out->rcont2[i],y_i1_value,y_i,dim,prec);
                
                //rcont3[i] = bspl = h * k1[i] - ydiff;
                _arb_vec_scalar_mul(v_s,k1,dim,h,prec);
                _arb_vec_sub(dense_out->rcont3[i],v_s,dense_out->rcont2[i],dim,prec);
                
                //rcont4[i] = ydiff - h*k4[i] - bspl;
                arb_add(s,x_i,h,prec); //fcn (n, xph, k5, pars, k4);
                func(k4, s, y_i1_value, dim, param, order, prec); //对 k4 进行更新，相当于下一次迭代的 k1
                
                _arb_vec_scalar_mul(v_s,k4,dim,h,prec);
                _arb_vec_sub(v_s,dense_out->rcont2[i],v_s,dim,prec);
                _arb_vec_sub(dense_out->rcont4[i],v_s,dense_out->rcont3[i],dim,prec);
                
                //rcont5[i] = d41*k1[i] + d46*k6[i] + d47*k7[i] + d48*k8[i] + 
                //            d49*k9[i] + d410*k10[i] + d411*k2[i] + d412*k3[i];
                _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_d401, prec);
                _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_d406, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_d407, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_d408, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k9, dim, ODEs_DOP853_coe_d409, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k10, dim, ODEs_DOP853_coe_d410, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k11, dim, ODEs_DOP853_coe_d411, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k12, dim, ODEs_DOP853_coe_d412, prec);
                _arb_vec_add(dense_out->rcont5[i],v_s,v_w,dim,prec);
                
                //rcont6[i] = d51*k1[i] + d56*k6[i] + d57*k7[i] + d58*k8[i] +
                //            d59*k9[i] + d510*k10[i] + d511*k2[i] + d512*k3[i];
                _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_d501, prec);
                _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_d506, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_d507, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_d508, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k9, dim, ODEs_DOP853_coe_d509, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k10, dim, ODEs_DOP853_coe_d510, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k11, dim, ODEs_DOP853_coe_d511, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k12, dim, ODEs_DOP853_coe_d512, prec);
                _arb_vec_add(dense_out->rcont6[i],v_s,v_w,dim,prec);
                
                
                //rcont7[i] = d61*k1[i] + d66*k6[i] + d67*k7[i] + d68*k8[i] +
                //            d69*k9[i] + d610*k10[i] + d611*k2[i] + d612*k3[i];
                _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_d601, prec);
                _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_d606, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_d607, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_d608, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k9, dim, ODEs_DOP853_coe_d609, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k10, dim, ODEs_DOP853_coe_d610, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k11, dim, ODEs_DOP853_coe_d611, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k12, dim, ODEs_DOP853_coe_d612, prec);
                _arb_vec_add(dense_out->rcont7[i],v_s,v_w,dim,prec);
                
                
                //rcont8[i] = d71*k1[i] + d76*k6[i] + d77*k7[i] + d78*k8[i] +
                //            d79*k9[i] + d710*k10[i] + d711*k2[i] + d712*k3[i];
                _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_d701, prec);
                _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_d706, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_d707, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_d708, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k9, dim, ODEs_DOP853_coe_d709, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k10, dim, ODEs_DOP853_coe_d710, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k11, dim, ODEs_DOP853_coe_d711, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k12, dim, ODEs_DOP853_coe_d712, prec);
                _arb_vec_add(dense_out->rcont8[i],v_s,v_w,dim,prec);
                
                
                /* the next three function evaluations */
                //k14 (k10)
                arb_mul(x_tem,ODEs_DOP853_coe_c14,h,prec);
                arb_add(x_tem,x_tem,x_i,prec);
                
                //y[i] + h * (a141*k1[i] + a147*k7[i] + a148*k8[i] + a149*k9[i] + a1410*k10[i] + a1411*k2[i] + 
                //            a1412*k3[i] + a1413*k4[i])
                _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_a1401, prec);
                _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_a1407, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_a1408, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k9, dim, ODEs_DOP853_coe_a1409, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k10, dim, ODEs_DOP853_coe_a1410, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k11, dim, ODEs_DOP853_coe_a1411, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k12, dim, ODEs_DOP853_coe_a1412, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k4, dim, ODEs_DOP853_coe_a1413, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                
                _arb_vec_scalar_mul(v_s, v_s, dim, h, prec);
                _arb_vec_add(y_tem, y_i, v_s, dim, prec);
                
                func(k14, x_tem, y_tem, dim, param, order, prec);
                
                //k15 (k2)
                arb_mul(x_tem,ODEs_DOP853_coe_c15,h,prec);
                arb_add(x_tem,x_tem,x_i,prec);
                
                //y[i] + h * (a151*k1[i] + a156*k6[i] + a157*k7[i] + a158*k8[i] +
                //            a1511*k2[i] + a1512*k3[i] + a1513*k4[i] + a1514*k10[i]);
                _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_a1501, prec);
                _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_a1506, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_a1507, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_a1508, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k11, dim, ODEs_DOP853_coe_a1511, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k12, dim, ODEs_DOP853_coe_a1512, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k4, dim, ODEs_DOP853_coe_a1513, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k14, dim, ODEs_DOP853_coe_a1514, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                
                _arb_vec_scalar_mul(v_s, v_s, dim, h, prec);
                _arb_vec_add(y_tem, y_i, v_s, dim, prec);
                
                func(k15, x_tem, y_tem, dim, param, order, prec);
                
                //k16 (k3)
                arb_mul(x_tem,ODEs_DOP853_coe_c16,h,prec);
                arb_add(x_tem,x_tem,x_i,prec);
                
                //y[i] + h * (a161*k1[i] + a166*k6[i] + a167*k7[i] + a168*k8[i] +
                //            a169*k9[i] + a1613*k4[i] + a1614*k10[i] + a1615*k2[i]);
                _arb_vec_scalar_mul(v_s, k1, dim, ODEs_DOP853_coe_a1601, prec);
                _arb_vec_scalar_mul(v_w, k6, dim, ODEs_DOP853_coe_a1606, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k7, dim, ODEs_DOP853_coe_a1607, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k8, dim, ODEs_DOP853_coe_a1608, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k9, dim, ODEs_DOP853_coe_a1609, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k4, dim, ODEs_DOP853_coe_a1613, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k14, dim, ODEs_DOP853_coe_a1614, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k15, dim, ODEs_DOP853_coe_a1615, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                
                _arb_vec_scalar_mul(v_s, v_s, dim, h, prec);
                _arb_vec_add(y_tem, y_i, v_s, dim, prec);
                
                func(k16, x_tem, y_tem, dim, param, order, prec);
                
                //利用 k14 k15 k16 再求 dense output 系数
                
                //rcont5[i] = h * (rcont5[i] + d413*k4[i] + d414*k10[i] + d415*k2[i] + d416*k3[i]);
                _arb_vec_scalar_mul(v_s, k4, dim, ODEs_DOP853_coe_d413, prec);
                _arb_vec_scalar_mul(v_w, k14, dim, ODEs_DOP853_coe_d414, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k15, dim, ODEs_DOP853_coe_d415, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k16, dim, ODEs_DOP853_coe_d416, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                
                _arb_vec_add(v_s,v_s,dense_out->rcont5[i],dim,prec);
                _arb_vec_scalar_mul(dense_out->rcont5[i], v_s, dim, h, prec);
                
                //rcont6[i] = h * (rcont6[i] + d513*k4[i] + d514*k10[i] + d515*k2[i] + d516*k3[i]);
                _arb_vec_scalar_mul(v_s, k4, dim, ODEs_DOP853_coe_d513, prec);
                _arb_vec_scalar_mul(v_w, k14, dim, ODEs_DOP853_coe_d514, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k15, dim, ODEs_DOP853_coe_d515, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k16, dim, ODEs_DOP853_coe_d516, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                
                _arb_vec_add(v_s,v_s,dense_out->rcont6[i],dim,prec);
                _arb_vec_scalar_mul(dense_out->rcont6[i], v_s, dim, h, prec);
                
                //rcont7[i] = h * (rcont7[i] + d613*k4[i] + d614*k10[i] + d615*k2[i] + d616*k3[i]);
                _arb_vec_scalar_mul(v_s, k4, dim, ODEs_DOP853_coe_d613, prec);
                _arb_vec_scalar_mul(v_w, k14, dim, ODEs_DOP853_coe_d614, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k15, dim, ODEs_DOP853_coe_d615, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k16, dim, ODEs_DOP853_coe_d616, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                
                _arb_vec_add(v_s,v_s,dense_out->rcont7[i],dim,prec);
                _arb_vec_scalar_mul(dense_out->rcont7[i], v_s, dim, h, prec);
                
                //rcont8[i] = h * (rcont8[i] + d713*k4[i] + d714*k10[i] + d715*k2[i] + d716*k3[i]);
                _arb_vec_scalar_mul(v_s, k4, dim, ODEs_DOP853_coe_d713, prec);
                _arb_vec_scalar_mul(v_w, k14, dim, ODEs_DOP853_coe_d714, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k15, dim, ODEs_DOP853_coe_d715, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                _arb_vec_scalar_mul(v_w, k16, dim, ODEs_DOP853_coe_d716, prec);
                _arb_vec_add(v_s,v_s,v_w,dim,prec);
                
                _arb_vec_add(v_s,v_s,dense_out->rcont8[i],dim,prec);
                _arb_vec_scalar_mul(dense_out->rcont8[i], v_s, dim, h, prec);
                
                
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
    arb_clear(sci);
    arb_clear(sum);
    arb_clear(sum_2);
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
    _arb_vec_clear(y_i1_diff_3,dim);
    _arb_vec_clear(y_i1_diff_5,dim);
    _arb_vec_clear(y_tem,dim);
    
    _arb_vec_clear(k1,dim);
    _arb_vec_clear(k2,dim);
    _arb_vec_clear(k3,dim);
    _arb_vec_clear(k4,dim);
    _arb_vec_clear(k5,dim);
    _arb_vec_clear(k6,dim);
    _arb_vec_clear(k7,dim);
    _arb_vec_clear(k8,dim);
    _arb_vec_clear(k9,dim);
    _arb_vec_clear(k10,dim);
    _arb_vec_clear(k11,dim);
    _arb_vec_clear(k12,dim);
    _arb_vec_clear(k14,dim);
    _arb_vec_clear(k15,dim);
    _arb_vec_clear(k16,dim);
}

