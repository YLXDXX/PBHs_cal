#include "header/06_draw_pic.h"

int C_r_prime_I(arb_t res, arb_t r, slong prec);
int C_r_prime_II(arb_t res, arb_t r, slong prec);


//画图
void draw_pic(char* comd_argv[], slong prec) // comd_argv 为命令行传递参数，可传递多个
{
    arb_t t,s,x_a,x_b,y_a,y_b,x,y,out_point,out_point_2;
    
    arb_init(t);
    arb_init(s);
    arb_init(x_a);
    arb_init(x_b);
    arb_init(y_a);
    arb_init(y_b);
    arb_init(x);
    arb_init(y);
    arb_init(out_point);
    arb_init(out_point_2);
    
    
    arb_set_str(x_a,"1E-10",prec); //一层循环用
    arb_set_str(x_b,"1E-6",prec);
    
    arb_set_str(y_a,"-0.6",prec); //二层循环用
    arb_set_str(y_b,"0.6",prec);
    
    slong out_number,number;
    
    number=5E2; //输出点的个数
    
    arb_ptr v_x_i,v_y_i;
    v_x_i=_arb_vec_init(number);
    v_y_i=_arb_vec_init(number);
    
    Get_interval_logspace_point(v_x_i,x_a,x_b,number,prec); //获取对数图上的间隔点
    //Get_interval_linspace_point(v_x_i,x_a,x_b,number,prec); //获取线性的间隔点
    
    //Get_interval_logspace_point(v_y_i,y_a,y_b,number,prec); //获取对数图上的间隔点
    Get_interval_linspace_point(v_y_i,y_a,y_b,number,prec); //获取线性的间隔点
    
    
    out_number=18; //输出数字的有效位数
    
    //写入文件
    FILE * fp;
    fp = fopen(Out_picture_file, "w"); //打开文件，a追加，w重新写入
    if( fp == NULL ) { //对文件打开操作进行判断
        printf("\n\nOpen Error: %s\t\n",Out_picture_file);perror("file");printf("\n");
        exit(-1);
    }
    
    //arb_set_str(PT_mu, "0.60", prec); //修改 PT_mu
    
    
    //一层输出
    for (slong i=0; i < number; i++)
    {
        arb_set(x,v_x_i+i);
        //arb_sub_ui(t,Upward_step_spectra_k_c,1,prec);
        //arb_log(t,x,prec);
        
        //power_spectrum(out_point,t,prec);
        //power_spectrum_non_Gaussian_f_Nl(out_point,x,prec);
        //zeta_Gauss_profile_n(out_point, x, 0, prec);
        //zeta_profile_n(out_point,x,0, prec); //ζ(r)
        //C_r_profile_n(out_point, x, 0, prec); //C(r)
        //Area_R(out_point, x, prec); //面积半径R(r)
        
        //C_r_prime_I(out_point,x,prec);
        //C_r_prime_II(out_point,x,prec);
        
        //PS相关
        //Abundance_beta_m_simplify(out_point,x,prec);
        //PS_abundance_beta_m(out_point, x, prec);
        //PS_abundance_f_m(out_point,x,prec);
        //Probability_zeta(out_point,x,NULL,0,prec);
        //Probability_C_l(out_point,x,prec);
        //Probability_C(out_point,x,prec);
        
        //考虑所有k模式，用δ谱计算连续谱
        //PS_abundance_beta_delta_k(out_point, x, prec); //计算某个k的β，临界坍缩的贡献都归到该k模式，传递值为ln(k)
        //PS_abundance_beta_delta_k_M(out_point,x,prec); //计算某个质量M(k)的β，考虑各个k的临界坍缩，传递值为ln(k)
        
        //诱导引力波
        //arb_exp(t,x,prec); //取指数后再传入
        //GW_power_spectra(out_point,eta,k,prec); //这里传入的k值未取对数
        //GW_current_energy_density(out_point,t,prec);
        //arb_set_str(x,"1",prec);
        //Func_GW_f_to_k(t, x, prec);//f nHz --> K Mpc^-1 
        //GW_current_energy_density_cuba(out_point,t,0,prec); //这里传入的k值未取对数
        //GW_current_energy_density_cuba(out_point_2,t,2,prec); //这里传入的k值未取对数
        //arb_printn(out_point, 50,0);printf("\n");
        //arb_set_str(x,"32",prec);
        //Func_GW_f_to_k(t, x, prec);//f nHz --> K Mpc^-1 
        //GW_current_energy_density_cuba(out_point,t,2,prec); //这里传入的k值未取对数
        //arb_printn(out_point, 50,0);printf("\n");
        //exit(0);
        
        //相对论自由度数
        //arb_exp(t,x,prec);
        //Effective_degrees_of_freedom_fit(out_point,out_point_2,t,"Gev",prec);
        
        
        arb_fprintn(fp,x,out_number,ARB_STR_NO_RADIUS); // 横坐标，对应参数值
        fprintf(fp, "\t");
        //arb_fprintn(fp,out_point_2,out_number,ARB_STR_NO_RADIUS); //out_point
        //fprintf(fp, "\t");
        arb_fprintn(fp,out_point,out_number,ARB_STR_NO_RADIUS); //out_point
        fprintf(fp, "\n");
        
        //printf("%ld/%ld\n",i,number);//进度显示
        print_progress(i,number); //进度条显示
    }
    
    
    /*
    //二层输出，例如，输出二维概率密度
    for (long int i=0; i < number; i++)
    {
        arb_set(x,v_x_i+i);
        
        for (long int j=0; j < number; j++)
        {
            arb_set(y,v_y_i+j);
            
            probability_gauss_2D(out_point,x,y,prec);
            
            arb_fprintn(fp,x,out_number,ARB_STR_NO_RADIUS); //out_point
            fprintf(fp, "\t");
            arb_fprintn(fp,y,out_number,ARB_STR_NO_RADIUS); //out_point
            fprintf(fp, "\t");
            arb_fprintn(fp,out_point,out_number,ARB_STR_NO_RADIUS); //out_point
            fprintf(fp, "\n");
            
            //printf("%ld/%ld\n",i*number+j,number*number);//进度显示
            print_progress(i*number+j,number*number); //进度条显示
        }
        
        //printf("%ld/%ld\n",i,number);//进度显示
    }
    */
    
    
    printf("\n输出完成\n");
    fclose(fp); //关闭文件
    
    arb_clear(t);
    arb_clear(s);
    arb_clear(x_a);
    arb_clear(x_b);
    arb_clear(y_a);
    arb_clear(y_b);
    arb_clear(x);
    arb_clear(y);
    arb_clear(out_point);
    arb_clear(out_point_2);
    
    _arb_vec_clear(v_x_i,number);
    _arb_vec_clear(v_y_i,number);
}


int C_r_prime_I(arb_t res, arb_t r, slong prec)
{
    arb_t t,s;
    
    arb_init(t);
    arb_init(s);
    
    // ζ’ + rζ’’ 模式
    zeta_profile_n(t,r,1,prec);
    
    zeta_profile_n(s,r,2,prec);
    
    arb_mul(s,s,r,prec);
    
    arb_add(res,t,s,prec);
    
    arb_clear(t);
    arb_clear(s);
    
    return 0;
}

int C_r_prime_II(arb_t res, arb_t r, slong prec)
{
    arb_t t,s;
    
    arb_init(t);
    arb_init(s);
    
    // 1 + rζ’ 
    zeta_profile_n(t,r,1,prec);
    
    arb_mul(t,t,r,prec);
    
    arb_add_ui(res,t,1,prec);
    
    arb_clear(t);
    arb_clear(s);
    
    return 0;
}



