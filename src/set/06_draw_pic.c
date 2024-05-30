#include "header/06_draw_pic.h"

int C_r_prime_I(arb_t res, arb_t r, slong prec);
int C_r_prime_II(arb_t res, arb_t r, slong prec);

//画图
void draw_pic(char* comd_argv, slong prec) // comd_argv 为命令行传递参数
{
    arb_t t,aa,bb,gap_x,gap_y,out_point,out_point_2,ay,by;
    
    arb_init(t);
    arb_init(aa); 
    arb_init(bb);
    arb_init(gap_x); //一层点间隔
    arb_init(gap_y); //二层点间隔
    arb_init(out_point);
    arb_init(out_point_2);
    
    arb_init(ay); 
    arb_init(by);
    
    arb_set_str(aa,"0.1",prec); //一层循环用
    arb_set_str(bb,"100",prec);
    
    /*
    arb_mul_ui(aa,Power_sigma, 30, prec); //aa=Ln_K_star-σ
    arb_sub(aa,Ln_K_star,aa,prec);
    arb_mul_ui(bb,Power_sigma, 500, prec);//bb=Ln_K_star+σ
    arb_add(bb,Ln_K_star,bb,prec);
    */
    
    arb_set_str(ay,"-1",prec); //二层循环用
    arb_set_str(by,"1",prec);
    
    slong out_number,number;
    
    number=1E3; //输出点的个数
    
    arb_sub(gap_x,bb,aa,prec); //x轴间隔
    arb_div_si(gap_x,gap_x,number,prec);
    
    arb_div_si(gap_x,gap_x,100,prec); //得到最小开始的点
    arb_add(aa,aa,gap_x,prec);
    arb_mul_si(gap_x,gap_x,100,prec);
    
    
    arb_sub(gap_y,by,ay,prec); //y轴间隔
    arb_div_si(gap_y,gap_y,number,prec);
    
    out_number=18; //输出数字的有效位数
    
    //写入文件
    FILE * fp;
    fp = fopen(Out_picture_file, "w"); //打开文件，a追加，w重新写入
    if( fp == NULL ) { //对文件打开操作进行判断
        printf("\n\nOpen Error: %s\t\n",Out_picture_file);perror("file");printf("\n");
        exit(-1);
    }
    
    //arb_set_str(Mu_2, "0.28", prec); //修改 Mu_2
    
    
    //一层输出
    for (long int i=1; i <= number; i++)
    {
        arb_fprintn(fp,aa,out_number,ARB_STR_NO_RADIUS); // 横坐标，对应参数值
        fprintf(fp, "\t");
        
        //zeta_Gauss_profile_n(out_point, aa, 0, prec);
        //zeta_profile_n(out_point,aa,0, prec); //ζ(r)
        //C_r_profile_n(out_point, aa, 0, prec); //C(r)
        
        //C_r_prime_I(out_point,aa,prec);
        //C_r_prime_II(out_point,aa,prec);
        
        //PS相关
        //Abundance_beta_m_simplify(out_point,aa,prec);
        //PS_abundance_beta_m(out_point, aa, prec);
        //PS_abundance_f_m(out_point,aa,prec);
        //Probability_zeta(out_point,aa,NULL,0,prec);
        //Probability_C_l(out_point,aa,prec);
        //Probability_C(out_point,aa,prec);
        
        //考虑所有k模式，用δ谱计算连续谱
        //PS_abundance_beta_delta_k(out_point, aa, prec); //计算某个k的β，临界坍缩的贡献都归到该k模式，传递值为ln(k)
        
        //诱导引力波
        //arb_exp(t,aa,prec); //取指数后再传入
        //GW_power_spectra(out_point,eta,k,prec); //这里传入的k值未取对数
        //GW_current_energy_density(out_point,t,prec);
        Func_GW_f_to_k(t, aa, prec);//f nHz --> K Mpc^-1 
        GW_current_energy_density_Omega_dim_2(out_point,t,prec); //这里传入的k值未取对数
        
        //相对论自由度数
        //arb_exp(t,aa,prec);
        //Effective_degrees_of_freedom_fit(out_point,out_point_2,t,"Gev",prec);
        
        
        //arb_fprintn(fp,out_point,out_number,ARB_STR_NO_RADIUS); //out_point
        //fprintf(fp, "\t");
        
        //PS_abundance_beta_delta_k_M(out_point,aa,prec); //计算某个质量M(k)的β，考虑各个k的临界坍缩，传递值为ln(k)
        
        //arb_printn(out_point, 50,0);printf("\n");
        
        //arb_fprintn(fp,out_point_2,out_number,ARB_STR_NO_RADIUS); //out_point
        //fprintf(fp, "\t");
        arb_fprintn(fp,out_point,out_number,ARB_STR_NO_RADIUS); //out_point
        fprintf(fp, "\n");
        
        arb_add(aa,aa,gap_x,prec); //更新点
        
        //printf("%ld/%ld\n",i,number);//进度显示
        print_progress(i,number); //进度条显示
    }
    
    
    /*
    //二层输出，例如，输出二维概率密度
    for (long int i=1; i <= number; i++)
    {
        arb_set_str(ay,"-1",prec); //每次循环，初始化第二层的值
        
        arb_div_si(gap_y,gap_y,100,prec);//得到第二层最小开始的点
        arb_add(ay,ay,gap_y,prec);
        arb_mul_si(gap_y,gap_y,100,prec);
        
        for (long int j=1; j <= number; j++)
        {
            probability_gauss_2D(out_point,aa,ay,prec);
            
            arb_fprintn(fp,aa,out_number,ARB_STR_NO_RADIUS); //out_point
            fprintf(fp, "\t");
            arb_fprintn(fp,ay,out_number,ARB_STR_NO_RADIUS); //out_point
            fprintf(fp, "\t");
            arb_fprintn(fp,out_point,out_number,ARB_STR_NO_RADIUS); //out_point
            fprintf(fp, "\n");
            
            arb_add(ay,ay,gap_y,prec); //更新内层点
            
            //printf("%ld/%ld\n",(i-1)*number+j,number*number);//进度显示
            print_progress((i-1)*number+j,number*number); //进度条显示
        }
        
        arb_add(aa,aa,gap_x,prec); //更新外层点
        
        //printf("%ld/%ld\n",i,number);//进度显示
    }
    */
    
    
    printf("\n输出完成\n");
    fclose(fp); //关闭文件
    
    arb_clear(t);
    arb_clear(aa);
    arb_clear(bb);
    arb_clear(gap_x);
    arb_clear(gap_y);
    arb_clear(out_point);
    arb_clear(out_point_2);
    arb_clear(ay);
    arb_clear(by);
    
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



