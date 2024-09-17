#include "routine_test.h"
#include <stdlib.h>
#include <arb_poly.h>
#include "fmpz.h"
#include <sys/time.h> 


//测试程序
void routine_test(slong prec)
{
    //运行时间
    //struct timeval start, end; 
    //long double time_taken; 
    
    arb_t s,t,w,x_start,x_end,error_abs,error_rel;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(x_start);
    arb_init(x_end);
    arb_init(error_abs);
    arb_init(error_rel);
    
    int dim=4;
    
    arb_ptr v_s,v_t,v_w,y_start;
    v_s=_arb_vec_init(dim);
    v_t=_arb_vec_init(dim);
    y_start=_arb_vec_init(dim);
    
    slong n=dim;
    v_w=_arb_vec_init(n);
    
    /*
    //函数插值测试
    
    gettimeofday(&start, NULL);
    for(slong i=0;i<n;i++)
    {
        arb_set_d(s,0.01*i+1.1); //x*sin(x)+cos(x)*e^(x)/(x^2-1)
        arb_sin(v_w+i,s,prec);
        arb_mul(v_w+i,v_w+i,s,prec);
        
        arb_sqr(t,s,prec);
        arb_sub_ui(t,t,1,prec);
        arb_exp(w,s,prec);
        arb_div(w,w,t,prec);
        arb_cos(t,s,prec);
        arb_mul(w,w,t,prec);
        arb_add(v_w+i,v_w+i,w,prec);
        
        arb_set(v_s+i,s);
    }
    
    gettimeofday(&end, NULL); 
    time_taken = (end.tv_sec  - start.tv_sec)  + (end.tv_usec  - start.tv_usec)  / 1e6; 
    printf("\nTime taken: %Lf seconds\n\n", time_taken); 
    
    
    arb_set_str(s,"3.13",prec);
    
    arb_poly_t p_s;
    arb_poly_init(p_s);
    
    // dim 个点全拟合，测试用
    //arb_poly_interpolate_newton(p_s,v_s,v_w,dim,prec); 
    //arb_poly_evaluate(w,p_s,s,prec);
    //arb_printn(w, 20,0);printf("\n");
    
    //分段拟合，实际用
    Interp_coe_t c_s; //用来存系数
    c_s=Interpolation_coe_init(dim);
    
    Interpolation_fit_func(w,s,v_s,v_w,c_s,dim,prec);
    arb_printn(w, 20,0);printf("\n");
    
    gettimeofday(&start, NULL);
    
    arb_set_str(s,"4.005",prec);
    
    Interpolation_fit_func(w,s,v_s,v_w,c_s,dim,prec);
    arb_printn(w, 20,0);printf("\n");
    //arb_poly_printd(p_s,5);
    
    arb_cos(w,s,prec); //x*sin(x)+cos(x)*e^(x)/(x^2-1)
    arb_exp(t,s,prec);
    arb_mul(w,w,t,prec);
    arb_sqr(t,s,prec);
    arb_sub_ui(t,t,1,prec);
    arb_div(w,w,t,prec);

    arb_sin(t,s,prec); 
    arb_mul(t,t,s,prec);
    arb_add(w,w,t,prec);
    
    arb_printn(w, 20,0);printf("\n");
    
    gettimeofday(&end, NULL); 
    time_taken = (end.tv_sec  - start.tv_sec)  + (end.tv_usec  - start.tv_usec)  / 1e6; 
    printf("\nTime taken: %Lf seconds\n\n", time_taken); 
    
    Interpolation_coe_clear(c_s,dim);
    
    exit(0);
    */
    
    // 初始条件
    /*
    phi_dot_0 = 1.0
    phi_0 = 1.0
    tau_0 = 0.0
    N_0 = 0.0
    H_0 = np.sqrt((1/6) * phi_dot_0**2 + V(phi_0) / 3)
    y0 = [phi_dot_0, phi_0, tau_0, N_0, H_0]
    */
    
    //微分方程测试
    
    arb_set_str(x_start,"0",prec); //初始条件
    arb_set_str(y_start,"0",prec);
    arb_set_str(y_start+1,"6.5",prec);
    arb_set_str(y_start+2,"0.0",prec);
    arb_set_str(y_start+3,"0.0",prec);
    
    
    arb_set_str(error_abs,"1E-15",prec);
    arb_set_str(error_rel,"1E-13",prec);
    
    arb_set_str(x_end,"1E7",prec); //需要足够大，方便后面计算
    
    
    /*
    ODEs_RFK45(v_s, Func_coupled_odes, dim, NULL, 0, //常微分方程组函数
               x_start, y_start, //给定初始条件
               x_end, //求出点 x_end 对应的函数值
               50, error, //num为迭代区间为[x_start,x_end]，将其分为几等份，从而给出初始步长
               prec);
    arb_printn(v_s, 50,0);printf("\n");
    arb_printn(v_s+1, 50,0);printf("\n");
    arb_printn(v_s+2, 50,0);printf("\n");
    arb_printn(v_s+3, 50,0);printf("\n");
    */
    
    slong num=1E5; //可能迭代数目的估计值
    
    G_Interp_coe_0=Interpolation_coe_init(num);
    G_Interp_coe_1=Interpolation_coe_init(num);
    G_Interp_coe_2=Interpolation_coe_init(num);
    G_Interp_coe_3=Interpolation_coe_init(num);
    G_Interp_coe_4=Interpolation_coe_init(num);
    
    
    
    
    //背景方程求解
    ODEs_DOPRI54_dense_t d_out;
    d_out=ODEs_DOPRI54_dense_init(num,dim);
    
    ODEs_DOPRI54(v_w,Func_background_coupled_odes,dim,NULL,0, //常微分方程组函数
                                     x_start, y_start, //给定初始条件
                                     x_end, //给定区间 [x_start, x_end]
                                     error_abs, error_rel, d_out,  //N为输出点个数，区间 N-1 等分，误差为相对精度
                                     prec);
    
    ///arb_printn(v_w+0, 50,0);printf("\n");
    arb_printn(v_w+0, 50,0);printf("\n");
    arb_printn(v_w+1, 50,0);printf("\n");
    arb_printn(v_w+2, 50,0);printf("\n");
    arb_printn(v_w+3, 50,0);printf("\n\n");
    
    if(0)
    {
        arb_printn(v_t+0, 50,0);printf("\n");
        arb_printn(v_w+0, 50,0);printf("\n");
    }
    
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
        //Interpolation_fit_func_odes_DOPRI54(out_y+j, out_x+j, d_out, 1, prec); //ϕ
        //Interpolation_fit_func_odes_DOPRI54(out_y+j, out_x+j, d_out, 0, prec); //ϕ'
        //Interpolation_fit_func_odes_DOPRI54(out_y+j, out_x+j, d_out, 3, prec); //N
        //Interpolation_fit_func_odes_DOPRI54(out_y+j, out_x+j, d_out, 2, prec); //τ
        
        
        //背景解 phi = phi_interp(t)
        Interpolation_fit_func_odes_DOPRI54(s, out_x+j, d_out, 1, prec);
        //背景解 phi_dot = phi_dot_interp(t)
        Interpolation_fit_func_odes_DOPRI54(t, out_x+j, d_out, 0, prec);
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
    arb_t t_ini,a_ini,Qphi_Real_t0,Qphi_Real_dot_t0;
    arb_t Qphi_Imag_t0,Qphi_Imag_dot_t0;
    
    arb_init(t_ini);
    arb_init(a_ini);
    arb_init(Qphi_Real_t0);
    arb_init(Qphi_Real_dot_t0);
    arb_init(Qphi_Imag_t0);
    arb_init(Qphi_Imag_dot_t0);
    
    arb_init(G_fourier_k);
    arb_set_str(G_fourier_k,"1E-1",prec);
    arb_set_str(t_ini,"0",prec);
    arb_set_str(a_ini,"1",prec);
    
    
    //得到 fk 对应的进入视界时间
    Func_get_time_k_enter(t_ini, s, G_fourier_k, t_ini, x_end, //在区间 [t_a, t_b] 内找根
                          a_ini, d_out, // a_i 为尺度因子的初始值
                          prec);
    //arb_printn(t_ini, 50,0);exit(0);
    
    
    //Qphi_Real_t0 = 1./(a_ini*np.sqrt(2.*fk)) //(np.cos(fk))/()
    arb_mul_si(s,G_fourier_k,2,prec);
    arb_sqrt(s,s,prec);
    arb_mul(s,s,a_ini,prec);
    arb_inv(Qphi_Real_t0,s,prec);
    
    //Qphi_Real_dot_t0 = 0.
    arb_set_str(Qphi_Real_dot_t0, "0", prec);
    
    //Qphi_Imag_t0 = 0. //(np.cos(fk))/()
    arb_set_str(Qphi_Imag_t0, "0", prec);
    
    //Qphi_Imag_dot_t0 = np.sqrt(fk)/(a_ini**2 * np.sqrt(2.))
    arb_one(s);
    arb_mul_si(s,s,2,prec);
    arb_sqrt(s,s,prec);
    arb_sqr(t,a_ini,prec);
    arb_mul(t,t,s,prec);
    arb_sqrt(s,G_fourier_k,prec);
    arb_div(Qphi_Imag_dot_t0,s,t,prec);
    
    // y_start = { Qphi_Real_t0, Qphi_Real_dot_t0, Qphi_Imag_t0, Qphi_Imag_dot_t0 }
    arb_set(y_start,Qphi_Real_t0);
    arb_set(y_start+1,Qphi_Real_dot_t0);
    arb_set(y_start+2,Qphi_Imag_t0);
    arb_set(y_start+3,Qphi_Imag_dot_t0);
    
    
    arb_set_str(x_end,"5E6",prec);
    
    ODEs_DOPRI54_dense_t d_out_2;
    d_out_2=ODEs_DOPRI54_dense_init(num,dim);
    
    ODEs_DOPRI54(v_s, Func_perturbation_phi_odes, dim, d_out, 0, //常微分方程组函数
               t_ini, y_start, //给定初始条件
               x_end, //求出点 x_end 对应的函数值
               error_abs, error_rel, NULL,
               prec);
    arb_printn(v_s, 50,0);printf("\n");
    arb_printn(v_s+1, 50,0);printf("\n");
    arb_printn(v_s+2, 50,0);printf("\n");
    arb_printn(v_s+3, 50,0);printf("\n\n");
    
    //exit(0);
    
    
    //利用扰动的解给出功率谱
    
    arb_t Qphi_Real_end,Qphi_Imag_end,H_end,phi_dot_end;
    arb_t N_end,phi_end,V_end,cos_fka,sin_fka;
    arb_init(Qphi_Real_end);
    arb_init(Qphi_Imag_end);
    arb_init(H_end);
    arb_init(phi_dot_end);
    arb_init(N_end);
    arb_init(phi_end);
    arb_init(V_end);
    arb_init(cos_fka);
    arb_init(sin_fka);
    
    slong nn=30;
    arb_ptr fi_k,fi_P;
    fi_k=_arb_vec_init(nn);
    fi_P=_arb_vec_init(nn);
    
    arb_set_str(s,"1E9",prec);
    arb_set_str(t,"1E10",prec);
    
    Get_interval_logspace_point(fi_k, s, t, nn, prec);
    
    arb_set_str(x_end,"5E6",prec);
    arb_zero(t_ini);
    arb_one(a_ini);
    
    arb_set_str(error_abs,"1E-30",prec);
    arb_set_str(error_rel,"1E-30",prec);
    
    for(slong i=0; i< nn; i++)
    {
        arb_zero(t_ini);
        arb_one(a_ini);
        
        arb_set(G_fourier_k,fi_k+i);
        
        //得到 fk 对应的进入视界时间，改变初始时间
        
        Func_get_time_k_enter(t_ini, s, fi_k+i, t_ini, x_end, //在区间 [t_a, t_b] 内找根
                              a_ini, d_out, // a_i 为尺度因子的初始值
                              prec);
        
        arb_exp(a_ini,s,prec);
        
        
        //cos(f_k/a_i*t_i)
        arb_div(s,fi_k+i,a_ini,prec);
        arb_mul(s,s,t_ini,prec);
        arb_cos(cos_fka,s,prec);
        arb_sin(sin_fka,s,prec);//sin(f_k/a_i*t_i)
        
        // Initial conditions for perturbation
        //arb_one(a_ini); //a_ini = 1.0
        
        //Qphi_Real_t0 = cos(f_k/a_i*t_i) / (a_ini * np.sqrt(2. * fk))
        //Qphi_Imag_t0 = sin(f_k/a_i*t_i) / (a_ini * np.sqrt(2. * fk))
        arb_mul_si(s,fi_k+i,2,prec);
        arb_sqrt(s,s,prec);
        arb_mul(s,s,a_ini,prec);
        
        arb_div(Qphi_Real_t0,cos_fka,s,prec);
        arb_div(Qphi_Imag_t0,sin_fka,s,prec);
        
        //Qphi_Real_dot_t0 = -sin(f_k/a_i*t_i)*np.sqrt(fk) / (a_ini**2 * np.sqrt(2.))
        //Qphi_Imag_dot_t0 = cos(f_k/a_i*t_i)*np.sqrt(fk) / (a_ini**2 * np.sqrt(2.))
        arb_set_ui(s,2);
        arb_sqrt(s,s,prec);
        arb_sqr(t,a_ini,prec);
        arb_mul(t,t,s,prec);
        arb_sqrt(s,fi_k+i,prec);
        arb_div(s,s,t,prec);
        
        arb_mul(Qphi_Real_dot_t0,sin_fka,s,prec);
        arb_neg(Qphi_Real_dot_t0,Qphi_Real_dot_t0);
        arb_mul(Qphi_Imag_dot_t0,cos_fka,s,prec);
        
        // y_start = { Qphi_Real_t0, Qphi_Real_dot_t0, Qphi_Imag_t0, Qphi_Imag_dot_t0 }
        arb_set(y_start,Qphi_Real_t0);
        arb_set(y_start+1,Qphi_Real_dot_t0);
        arb_set(y_start+2,Qphi_Imag_t0);
        arb_set(y_start+3,Qphi_Imag_dot_t0);
        
        //arb_printn(t_ini, 50,0);printf("\n");exit(0);
        //arb_printn(G_fourier_k, 50,0);printf("\n\n");exit(0);
        //arb_printn(H_end, 50,0);printf("\n\n");
        
        ODEs_DOPRI54(v_s, Func_perturbation_phi_odes, dim, d_out, 0, //常微分方程组函数
                     t_ini, y_start, //给定初始条件
                     x_end, //求出点 x_end 对应的函数值
                     error_abs, error_rel, d_out_2,
                     prec);
        
        /*
        //输出实部和虚部图像
        slong out_num=5000;
        arb_ptr out_x,out_y;
        out_x=_arb_vec_init(out_num);
        out_y=_arb_vec_init(out_num);
        
        Get_interval_logspace_point(out_x, t_ini, x_end, out_num, prec);
        for(slong j=0; j< out_num;j++)
        {
            Interpolation_fit_func_odes_DOPRI54(out_y+j, out_x+j, d_out_2, 0, prec);
        }
        Vector_point_output_to_file(out_x, out_y, out_num, 'w'); //a追加，w重新写入
        
        arb_printn(out_y, 50,0);printf("\n");
        arb_printn(out_y+out_num-1, 50,0);printf("\n");
        exit(0);
        //输出实部和虚部图像 完
        */
        
        // Calculate P_R at the final time
        arb_set(Qphi_Real_end,v_s+0);
        arb_set(Qphi_Imag_end,v_s+2);
        
        
        
        //背景解 N = N_interp(t)
        Interpolation_fit_func_odes_DOPRI54(N_end, x_end, d_out, 3,  prec);
        
        //背景解 phi = phi_interp(t)
        Interpolation_fit_func_odes_DOPRI54(phi_end, x_end, d_out, 1, prec);
        
        //背景解 phi_dot = phi_dot_interp(t)
        Interpolation_fit_func_odes_DOPRI54(phi_dot_end, x_end, d_out, 0, prec);
        
        Func_V_phi(V_end,phi_end,prec); //V = V_phi(phi)
        
        //背景解 H = H_interp(t), 这里不用插值
        //H = np.sqrt((1./6.) * phi_dot**2 + V / 3.)
        arb_sqr(s,phi_dot_end,prec);
        arb_div_ui(s,s,6,prec);
        arb_div_ui(w,V_end,3,prec);
        arb_add(s,s,w,prec);
        arb_sqrt(H_end,s,prec); //得到 H
        
        //P_R = ((fk**3) / (2 * np.pi**2)) * ((H_end*Qphi_Real_end/phi_dot_end)**2 + (H_end*Qphi_Imag_end/phi_dot_end)**2)
        
        //((fk**3) / (2 * np.pi**2))
        arb_sqr(t,Pi,prec);
        arb_mul_si(t,t,2,prec);
        arb_pow_ui(w,fi_k+i,3,prec);
        arb_div(w,w,t,prec);
        
        
        //(H_end*Qphi_Real_end/phi_dot_end)**2
        arb_mul(s,H_end,Qphi_Real_end,prec);
        arb_div(s,s,phi_dot_end,prec);
        arb_sqr(s,s,prec);
        
        //(H_end*Qphi_Imag_end/phi_dot_end)**2
        arb_mul(t,H_end,Qphi_Imag_end,prec);
        arb_div(t,t,phi_dot_end,prec);
        arb_sqr(t,t,prec);
        
        arb_add(s,s,t,prec);
        arb_mul(fi_P+i,w,s,prec);
        
        printf("%li\n",i);
        //arb_printn(y_start+2, 50,0);printf("\n");
        arb_printn(fi_P+i, 50,0);printf("\n\n");
    }
    
    Vector_point_output_to_file(fi_k, fi_P, nn, 'w'); //a追加，w重新写入
    
    exit(0);
    
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
    /*
    // Func_test_ODEs_func_04 精确解为 y=(1+1/6*x^3)*e^{x}
    arb_pow_ui(s,x_end,3,prec);
    arb_div_ui(s,s,6,prec);
    arb_add_ui(s,s,1,prec);
    arb_exp(t,x_end,prec);
    arb_mul(s,s,t,prec);
    */
    
    //arb_printn(s, 50,0);printf("\n");
    
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
 
