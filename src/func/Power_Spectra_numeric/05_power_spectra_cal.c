#include "Inflation_odes_slove.h"
#include <stdlib.h> 
#include <omp.h>


//功率谱的数值计算
void Inflation_power_spectra_numeric_cal(slong prec)
{
    // dense output 的类型为 Inflation_dense_t, 初始化 Inflation_dense_init, 清理 Inflation_dense_clear
    // dense output 插值拟合 Inflation_interp_fit_func_odes
    // 扰动ODEs参数 Inflation_perturb_ODEs_param_t,
    // 扰动ODEs参数: 初始化 Inflation_perturb_ODEs_param_init, 清理 Inflation_perturb_ODEs_param_clear
    // 常微分方程求解 Inflation_ODEs_solver
    
    Inflation_set_model_parameters(prec); //计算前，先设定模型参数
    
    arb_t s,t,w,N,x_start,x_end,error_abs,error_rel;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(N);
    arb_init(x_start);
    arb_init(x_end);
    arb_init(error_abs);
    arb_init(error_rel);
    
    arb_t t_ini,a_ini,fk; //解扰动方程用
    arb_init(t_ini);
    arb_init(a_ini);
    arb_init(fk);
    
   
    
    slong dim=4; //方程组维数
    slong num=1E5; //可能迭代数目的估计值
    
    arb_ptr v_s,v_w,y_start;
    v_s=_arb_vec_init(dim);
    v_w=_arb_vec_init(dim);
    y_start=_arb_vec_init(dim);
    
    Inflation_dense_t d_out; //用于密度输出
    d_out=Inflation_dense_init(num,dim);
    
    arb_set_str(x_start,"0",prec); //初始条件, t
    arb_set_str(y_start,"0",prec); // dot(ϕ)
    arb_set_str(y_start+1,"6.5",prec); // ϕ
    arb_set_str(y_start+2,"0.0",prec); //τ
    arb_set_str(y_start+3,"0.0",prec); //N
    
    arb_set_str(error_abs,"1E-15",prec); //误差
    arb_set_str(error_rel,"1E-13",prec);
    
    arb_set_str(x_end,"1E7",prec); //需要足够大，方便后面计算
    
    
    //背景方程求解
    Inflation_ODEs_solver(v_w,Inflation_background_coupled_odes,dim,NULL,0, //常微分方程组函数
                          x_start, y_start, //给定初始条件
                          x_end, //给定区间 [x_start, x_end]
                          error_abs, error_rel, d_out,  //N为输出点个数，区间 N-1 等分，误差为相对精度
                          prec);
    
    
    arb_printn(v_w+0, 50,0);printf("\n");
    arb_printn(v_w+1, 50,0);printf("\n");
    arb_printn(v_w+2, 50,0);printf("\n");
    arb_printn(v_w+3, 50,0);printf("\n\n");
    
    
    /*
    //背景解作图
    slong out_num=1E3;
    arb_ptr out_x,out_y;
    out_x=_arb_vec_init(out_num);
    out_y=_arb_vec_init(out_num);
    
    arb_set_str(s,"1",prec);
    arb_set_str(t,"1E7",prec);
    Get_interval_logspace_point(out_x, s, t, out_num, prec);
    //Get_interval_linspace_point(out_x, s, t, out_num, prec);
    for(slong j=0; j< out_num;j++)
    {
        //Func_V_phi_p(out_y+j, , prec);
        //Inflation_interp_fit_func_odes(out_y+j, out_x+j, d_out, 1, prec); //ϕ
        //Inflation_interp_fit_func_odes(out_y+j, out_x+j, d_out, 0, prec); //ϕ'
        //Inflation_interp_fit_func_odes(out_y+j, out_x+j, d_out, 3, prec); //N
        //Inflation_interp_fit_func_odes(out_y+j, out_x+j, d_out, 2, prec); //τ
        
        
        //背景解 phi = phi_interp(t)
        Inflation_interp_fit_func_odes(s, out_x+j, d_out, 1, prec);
        //背景解 phi_dot = phi_dot_interp(t)
        Inflation_interp_fit_func_odes(t, out_x+j, d_out, 0, prec);
        Func_V_phi(w,s,prec); //V = V_phi(phi)
        
        //背景解 H = H_interp(t), 这里不用插值
        //H = np.sqrt((1./6.) * phi_dot**2 + V / 3.)
        arb_sqr(s,t,prec);
        arb_div_ui(s,s,6,prec);
        arb_div_ui(w,w,3,prec);
        arb_add(s,s,w,prec);
        arb_sqrt(out_y+j,s,prec); //得到 H
        
    }
    Vector_point_output_to_file(out_x, out_y, out_num, 'w'); //a追加，w重新写入
    exit(0);
    //背景解作图 完
    */
    
    
    //扰动方程求解
    arb_set_str(fk,"1E-1",prec); //初始条件
    arb_set_str(t_ini,"0",prec);
    arb_set_str(a_ini,"1",prec);
    
    
    //得到 fk 对应的进入视界时间, 输出对应的时间 和 N 值
    //改变初始时间，初始尺度因子
    Inflation_get_time_k_enter(t_ini, N, fk, t_ini, x_end, //在区间 [t_a, t_b] 内找根
                          a_ini, d_out, // a_i 为尺度因子的初始值
                          prec);
    arb_exp(a_ini,N,prec);
    //arb_printn(t_ini, 50,0);exit(0);
    
    
    //求扰动方程的初始条件
    Inflation_get_perturb_oeds_initial_condition(y_start, t_ini, a_ini, fk, prec);
    
    arb_set_str(x_end,"5E6",prec);
    
    //arb_set_str(error_abs,"1E-23",prec);
    //arb_set_str(error_rel,"1E-23",prec);
    //参数传入
    Inflation_perturb_ODEs_param_t param_p;
    param_p=Inflation_perturb_ODEs_param_init(d_out,fk); //背景的dense output 和 当前的 fk
    
    Inflation_ODEs_solver(v_s, Inflation_perturbation_phi_odes, dim, param_p, 0, //常微分方程组函数
                          t_ini, y_start, //给定初始条件
                          x_end, //求出点 x_end 对应的函数值
                          error_abs, error_rel, NULL,
                          prec);
    
    arb_printn(v_s, 50,0);printf("\n");
    arb_printn(v_s+1, 50,0);printf("\n");
    arb_printn(v_s+2, 50,0);printf("\n");
    arb_printn(v_s+3, 50,0);printf("\n\n");
    
    //exit(0);
    
    slong nn=100;
    arb_ptr fi_k,fi_P;
    fi_k=_arb_vec_init(nn);
    fi_P=_arb_vec_init(nn);
    
    arb_set_str(s,"1E2",prec);
    arb_set_str(t,"3E5",prec);
    
    Get_interval_logspace_point(fi_k, s, t, nn, prec);
    
    arb_set_str(x_end,"5E6",prec);
    arb_zero(t_ini);
    arb_one(a_ini);
    
    arb_set_str(error_abs,"1E-25",prec);
    arb_set_str(error_rel,"1E-25",prec);
    
    
    //改造for循环后，可以使用 多线程加速，采用 OpenMP 方法
    #pragma omp parallel for num_threads(10)
    for(slong i=0; i< nn; i++)
    {
        //多线程并发，需保证每次循环中间使用各变量独立
        arb_t tt_ini,aa_ini,NN;
        arb_init(tt_ini);
        arb_init(aa_ini);
        arb_init(NN);
        
        arb_ptr yy_start,yy_end;
        yy_start=_arb_vec_init(dim);
        yy_end=_arb_vec_init(dim);
        
        arb_zero(tt_ini);//初始值
        arb_one(aa_ini);
        
        Inflation_perturb_ODEs_param_t param_pp; //参数传递
        param_pp=Inflation_perturb_ODEs_param_init(d_out, fi_k+i);
        
        
        //得到 fk 对应的进入视界时间，改变初始时间，初始尺度因子
        Inflation_get_time_k_enter(tt_ini, NN, fi_k+i, tt_ini, x_end, //在区间 [t_a, t_b] 内找根
                              aa_ini, d_out, // aa_i 为尺度因子的初始值
                              prec);
        arb_exp(aa_ini,NN,prec);
        
        
        //求扰动方程的初始条件
        Inflation_get_perturb_oeds_initial_condition(yy_start, tt_ini, aa_ini, fi_k+i, prec);
        
        //扰动方程求解
        Inflation_ODEs_solver(yy_end, Inflation_perturbation_phi_odes, dim, param_pp, 0, //常微分方程组函数
                              tt_ini, yy_start, //给定初始条件
                              x_end, //求出点 x_end 对应的函数值
                              error_abs, error_rel, NULL,
                              prec);
        
        /*
        //输出*实部和虚部图像
        slong out_num=5000;
        arb_ptr out_x,out_y;
        out_x=_arb_vec_init(out_num);
        out_y=_arb_vec_init(out_num);
        
        Get_interval_logspace_point(out_x, t_ini, x_end, out_num, prec);
        for(slong j=0; j< out_num;j++)
        {
            Inflation_interp_fit_func_odes(out_y+j, out_x+j, d_out_2, 0, prec);
        }
        Vector_point_output_to_file(out_x, out_y, out_num, 'w'); //a追加，w重新写入
        
        arb_printn(out_y, 50,0);printf("\n");
        arb_printn(out_y+out_num-1, 50,0);printf("\n");
        exit(0);
        //输出实部和虚部图像 完
        */
        
        //计算 fk 模式的功率谱
        Inflation_power_spectra_cal_at_fk(fi_P+i, fi_k+i, x_end, yy_end, d_out, prec);
        
        printf("%li\n",i);
        //arb_printn(y_start+2, 50,0);printf("\n");
        //arb_printn(fi_P+i, 50,0);printf("\n\n");
        
        
        arb_clear(tt_ini);
        arb_clear(aa_ini);
        arb_clear(NN);
        
        _arb_vec_clear(yy_start,dim);
        _arb_vec_clear(yy_end,dim);
        
        //清理释放结构体
        Inflation_perturb_ODEs_param_clear(param_pp);
    }
    
    //将功率谱的结果输出到文件
    Vector_point_output_to_file(fi_k, fi_P, nn, 'w'); //a追加，w重新写入
    exit(0);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(N);
    arb_clear(x_start);
    arb_clear(x_end);
    arb_clear(error_abs);
    arb_clear(error_rel);
    
    arb_clear(t_ini);
    arb_clear(a_ini);
    arb_clear(fk);
    
    _arb_vec_clear(v_s,dim);
    _arb_vec_clear(v_w,dim);
    _arb_vec_clear(y_start,dim);
    
    Inflation_dense_clear(d_out, num, dim);
    
    Inflation_perturb_ODEs_param_clear(param_p);
}
 
