#include "C_l_probability.h" 
#include <arb_mat.h> //矩阵
#include <stdlib.h>

//二维高斯概率分布P(X, Y)
int probability_gauss_2D(arb_t res, const arb_t x, const arb_t y, slong prec)
{
    //局部变量设定
    arb_t t,s,det;
    arb_mat_t V,T,S,Sigma;
    
    //初始化参数
    arb_init(t);
    arb_init(s);
    arb_init(det);
    
    //初始化矩阵
    arb_mat_init(V, 2, 1); // 二行一列
    arb_mat_init(T, 1, 2);
    arb_mat_init(S, 1, 1);
    arb_mat_init(Sigma, 2, 2);//二阶协方差矩阵
    
    arb_set(arb_mat_entry(V, 0, 0),x);
    arb_set(arb_mat_entry(V, 1, 0),y);
    
    //设定协方差阵的值
    arb_set(arb_mat_entry(Sigma, 0, 0),PS_Sigma_XX);
    arb_set(arb_mat_entry(Sigma, 0, 1),PS_Sigma_XY);
    arb_set(arb_mat_entry(Sigma, 1, 0),PS_Sigma_YX);
    arb_set(arb_mat_entry(Sigma, 1, 1),PS_Sigma_YY);
    
    //arb_mat_printd(Sigma, 3);
    
    arb_mat_transpose(T,V); //V矩阵转置
    arb_mat_det(det,Sigma,prec); //矩阵行列式
    arb_mat_inv(Sigma,Sigma,prec); //矩阵求逆
    
    //指数部分 V^T * Sigma^{-1} * V
    arb_mat_mul(T,T,Sigma,prec); //矩阵乘积
    arb_mat_mul(S,T,V,prec);
    arb_set(s,arb_mat_entry(S, 0, 0)); // S为1x1矩阵，转为实数
    arb_div_si(s,s,2,prec);
    arb_neg(s,s);
    arb_exp(s,s,prec);
    
    
    //前面系数部分 1/[2π*sqrt(det)]
    arb_abs(det,det); //开根号前绝对值
    arb_sqrt(t,det,prec);
    arb_mul(t,t,Pi_2,prec);
    arb_inv(t,t,prec); //Sets z to 1/𝑥
    
    arb_mul(res,s,t,prec);
    
    //完成计算，释放
    arb_clear(t);
    arb_clear(s);
    arb_clear(det);
    arb_mat_clear(V);
    arb_mat_clear(T);
    arb_mat_clear(S);
    arb_mat_clear(Sigma);
    
    return 0;
}



// C_ℓ=-4/3 * [ J_1(Y)*X + 2*J_2(Y)*(Σ_XY) ]
// C_ℓ 展开系数 J_1
int interior_C_l_J_1(arb_t res, const arb_t Y, slong prec)
{
    arb_t s;
    
    arb_init(s);
    
    //针对各种不同的 ζ 类型 分开讨论， Y=ζ_G
    switch(Zeta_type) 
    {
        case gaussian_type : 
            //对于纯高斯的情况 ζ=ζ_G --> J_1=1
            arb_one(res);
            break;
            
        case exponential_tail_type :
            //此时 ζ=-1/3 * ln(1-3*ζ_G)  一阶导 --> J_1=1/(1-3*ζ_G)
            Non_Gaussianity_exponential_tail_n(s,Y,1,prec); //一阶导，调用前面函数
            arb_set(res,s);
            break;
            
        case up_step_type :
            // up_step 非高斯性的函数为 R=-2/|h|*[sqrt(1-|h|R_G)-1]
            //即 y=-2/h*[sqrt(1-h*x)-1]
            Non_Gaussianity_up_step_n(s,Y,1,prec); //一阶导，调用前面函数
            arb_set(res,s);
            break;
            
        case power_expansion_type :
            //power_expansion 非高斯性的函数为 ζ(r)=ζ_G(r) + 3/5 * f_NL * ζ_G^2 + 9/25*g_NL*ζ_G^3
            Non_Gaussianity_power_expansion_n(s,Y,1,prec); //一阶导，调用前面函数
            arb_set(res,s);
            break;
            
        case narrow_step_1_type :
            //有限宽 step 非高斯性，扰动 1 阶
            Non_Gaussianity_narrow_1_up_step_n(s,Y,1,prec); //一阶导，调用前面函数
            arb_set(res,s);
            break;
            
        case narrow_step_1_2_type :
            //有限宽 step 非高斯性，扰动 1+2 阶
            Non_Gaussianity_narrow_1_2_up_step_n(s,Y,1,prec); //一阶导，调用前面函数
            arb_set(res,s);
            break;
            
        default :
            printf("PS_Method -> C_l_probability -> 在 interior_C_l_J_1 -> Zeta_type不正确\n" );
            exit(1);
    }
    
    //完成计算，释放
    arb_clear(s);
    
    return 0;
}

// C_ℓ 展开系数 J_2
int interior_C_l_J_2(arb_t res, const arb_t Y, slong prec)
{
    //针对各种不同的 ζ 类型 分开讨论， Y=ζ_G
    //当 ζ=f(ζ_G) 只是 ζ_G的函数 时， J_2=0
    //当前，对于各种 Zeta_type 类型，都设为零
    
    arb_zero(res); //设为零
    
    return 0;
}


//计算概率 P(C_l) 用 
int interior_probability_C_l(arb_t res, const arb_t Y, void* C_l, const slong order, slong prec)
{
    arb_t t,s,u,w;
    
    //初始化参数
    arb_init(t);
    arb_init(s);
    arb_init(u);
    arb_init(w);
    
    //其中 C_l 的值，通过参数的形式传入， Y=ζ_G
    
    //左边系数部分 3/[4*|(J_1(Y))|]
    interior_C_l_J_1(t,Y,prec);
    arb_inv(t,t,prec);
    arb_set(w,t); //留作右边备用
    arb_abs(t,w);//取绝对值
    arb_mul_si(t,t,3,prec);
    arb_div_si(t,t,4,prec);
    
    
    //右边二维概率分布部分
    arb_mul_si(s,C_l,3,prec);
    arb_div_si(s,s,4,prec);
    
    interior_C_l_J_2(u,Y,prec);
    arb_mul(u,u,PS_Sigma_XY,prec);
    arb_mul_si(u,u,2,prec);
    arb_add(u,u,s,prec);
    
    arb_mul(u,u,w,prec);
    arb_neg(u,u);
    
    arb_set(w,Y);
    
    probability_gauss_2D(s,u,w,prec);
    
    
    arb_mul(res,s,t,prec);
    
    //完成计算，释放
    arb_clear(t);
    arb_clear(s);
    arb_clear(u);
    arb_clear(w);
    return 0;
}

int interior_probability_C_l_Y_root(arb_t res, const arb_t Y, void* parameter, const slong order, slong prec)
{
    arb_t t,s;
    
    //初始化参数
    arb_init(t);
    arb_init(s);
    
    //求根方程： dζ/dζ_G*Y+A*C_l=0
    
    //需传入结构体 Find_root_delta_C_l_Y 来获取参数 A 和 C_l
    struct Find_root_delta_C_l_Y *Root_Y_parameter; //这里不需要手动分配，只需将其指向传入的指针即可
    Root_Y_parameter=parameter;
    
    
    interior_C_l_J_1(t,Y,prec); //dζ/dζ_G
    arb_mul(t,t,Y,prec); //dζ/dζ_G*Y
    
    arb_mul(s, Root_Y_parameter->A, Root_Y_parameter->C_l, prec); //A*C_l
    
    arb_add(res,t,s,prec);
    
    //完成计算，释放
    arb_clear(t);
    arb_clear(s);
    return 0;
}

//多个Y时，求每个Y对应的概率
void interior_probability_C_l_P_each_Y(arb_t res, const arb_t Y, const arb_t A, const arb_t C_l, slong prec )
{
    arb_t s,t,zeta_1,zeta_2;
    arb_init(s);
    arb_init(t);
    arb_init(zeta_1);
    arb_init(zeta_2);
    
    //对应概率： P(C_l)=P(Y)*|A|*|ζ''Y+ζ'|^(-1)  //取绝对值
    // or P(C_l)==P(Y)*|C_l/(Y*ζ')*(ζ''Y+ζ')|^(-1) //此式适用面更广
    
    //为求|C_l/(Y*ζ')*(ζ''Y+ζ')|^(-1)，针对各种不同的 ζ 类型 分开讨论， Y=ζ_G
    switch(Zeta_type) 
    {
        case exponential_tail_type :
            Non_Gaussianity_exponential_tail_n(zeta_1,Y,1,prec); //一阶导，调用前面函数
            Non_Gaussianity_exponential_tail_n(zeta_2,Y,2,prec); //二阶导，调用前面函数
            break;
            
        case up_step_type :
            Non_Gaussianity_up_step_n(zeta_1,Y,1,prec); //一阶导，调用前面函数
            Non_Gaussianity_up_step_n(zeta_2,Y,2,prec); //二阶导，调用前面函数
            break;
            
        case power_expansion_type :
            Non_Gaussianity_power_expansion_n(zeta_1,Y,1,prec); //一阶导，调用前面函数
            Non_Gaussianity_power_expansion_n(zeta_2,Y,2,prec); //二阶导，调用前面函数
            break;
            
        case narrow_step_1_type :
            Non_Gaussianity_narrow_1_up_step_n(zeta_1,Y,1,prec); //一阶导，调用前面函数
            Non_Gaussianity_narrow_1_up_step_n(zeta_2,Y,2,prec); //二阶导，调用前面函数
            break;
            
        case narrow_step_1_2_type :
            Non_Gaussianity_narrow_1_2_up_step_n(zeta_1,Y,1,prec); //一阶导，调用前面函数
            Non_Gaussianity_narrow_1_2_up_step_n(zeta_2,Y,2,prec); //二阶导，调用前面函数
            break; 
        default :
            printf("PS_Method -> C_l_probability -> interior_probability_C_l_P_each_Y -> Zeta_type不正确\n" );
            exit(1);
    }
    
    arb_mul(t,zeta_2,Y,prec); //|C_l/(Y*ζ')*(ζ''Y+ζ')|^(-1)
    arb_add(t,t,zeta_1,prec);
    //有一些模型，有取值范围的限制，超过取值范围后，ζ''和ζ'的取值均为零
    if(arb_is_zero(t)) //t=0
    {
        //printf("超出范围\n");
        arb_zero(res);
    }else
    {
        arb_mul(t,t,C_l,prec);
        arb_div(t,t,Y,prec);
        arb_div(t,t,zeta_1,prec);
        
        arb_inv(t,t,prec); //取倒数
        arb_abs(t,t); //取绝对值
        
        //P(Y)
        interior_probability_gauss(s,Y,prec);
        
        arb_mul(res,s,t,prec);
    }
    
    
    arb_clear(t);
    arb_clear(s);
    arb_clear(zeta_1);
    arb_clear(zeta_2);
}

//计算 C_ℓ 的概率密度分布 P(C_l)
int Probability_C_l(arb_t res, const arb_t C_l, slong prec)
{
    arb_t s,t,w,A,Y,t_cl,P_each_Y;
    
    //初始化参数
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(A);
    arb_init(Y);
    arb_init(t_cl);
    arb_init(P_each_Y);
    
    arb_set(t_cl,C_l);
    
    if( Zeta_type==gaussian_type )
    {
        //在高斯情况下，连续谱，可以不用求积分，P(C_l)的分布是高斯分布
        //δ谱也对应一个高斯分布，与连续谱的情况一样
        //可直接用高斯分布求解
        //Σ_{C_l}=(4/3)^2 * Σ_XX
        arb_one(s);
        arb_mul_ui(s,s,16,prec);
        arb_div_ui(s,s,9,prec);
        arb_mul(s,s,PS_Sigma_XX,prec);
        
        //高斯分布，指数部分
        arb_sqr(t,C_l,prec);
        arb_div_ui(t,t,2,prec);
        arb_div(t,t,s,prec);
        arb_neg(t,t);
        arb_exp(t,t,prec);
        
        //高斯分布，系数部分
        arb_mul(w,Pi_2,s,prec);
        arb_sqrt(w,w,prec);
        
        arb_div(res,t,w,prec);
        
    }else
    {
        switch(Power_spectrum_type)
        {
            case lognormal_type : //此几种情况可合并，在最后加 break 即可
            case power_law_type :
            case broken_power_law_type :
            case box_type :
            case link_cmb_type :
                //连续谱，非高斯情况，使用积分求解
                Integration_arb(w, interior_probability_C_l, t_cl, 0, 
                                PS_Int_P_C_l_min, PS_Int_P_C_l_max,PS_Int_P_C_l_precision,
                                Integration_iterate_min,Integration_iterate_max, prec);
                
                arb_set(res,w);
                
                break;
            case delta_type : 
                //利用附录中对于δ情况的推导求解
                //δ情况下，P(X,Y)中由于X和Y线性相关，退化为一元情况
                //需要通过 C_l 反解出 Y
                //有多个根的情况下，是将各个根对应的概率加起来
                
                //对于给定的C_l需要求出Y，可能有多个根
                //dζ/dζ_G*Y+3/4*C_l*[sin(x)/( xcos(x)-sin(x) )]=0
                // or dζ/dζ_G*Y+3/4*C_l*Σ_YY/Σ_XY=0 //此式适用面更广
                //设常数A=3/4*[sin(x)/( xcos(x)-sin(x) )]
                //dζ/dζ_G*Y+A*C_l=0
                
                arb_div(s,PS_Sigma_YY,PS_Sigma_XY,prec); //系数A
                arb_mul_ui(s,s,3,prec);
                arb_div_ui(A,s,4,prec);
                
                //反解Y，求根传参数
                //这里，对于结构体 Find_root_delta_C_l_Y 需手动分配内存
                struct Find_root_delta_C_l_Y *Root_Y_parameter = (struct Find_root_delta_C_l_Y *)calloc(1,sizeof(struct Find_root_delta_C_l_Y));
                
                arb_init(Root_Y_parameter->C_l);//使用arb_t变量前初始化
                arb_init(Root_Y_parameter->A);
                
                arb_set(Root_Y_parameter->C_l,t_cl); //设定求根参数
                arb_set(Root_Y_parameter->A,A);
                
                int root_num; //根的个数
                arb_ptr muil_r; //存储多个根
                arb_ptr* m_r; //改变muil_r指针指向的地址，需要一个指向该指针的指针
                m_r=&muil_r;
                
                //arb_printn(A, 50,0);printf("\n");
                //arb_printn(t_cl, 50,0);printf("\n");
                
                switch(Zeta_type)
                {
                    case exponential_tail_type : //此几种情况可合并，在最后加 break 即可
                    case up_step_type :
                    case power_expansion_type :
                    case narrow_step_1_type :
                    case narrow_step_1_2_type :
                        
                        //此时 P(C_l)=P(Y)*|A|*|ζ''Y+ζ'|^(-1)  //取绝对值
                        //需要反解出Y，这里用数值解法
                        //dζ/dζ_G*Y+A*C_l=0
                        
                        root_num=Find_interval_multi_root(m_r,interior_probability_C_l_Y_root,Root_Y_parameter, 0,
                                                          PS_Root_C_l_to_Y_min,PS_Root_C_l_to_Y_max,PS_Root_C_l_to_Y_precision,
                                                          PS_Root_C_l_to_Y_num,prec);
                        //printf("根的个数为: %i\n",root_num);
                        //将每个根对应的根率加起来
                        arb_zero(P_each_Y);
                        for(int root_i=0; root_i<root_num; root_i++)
                        {
                            interior_probability_C_l_P_each_Y(w,muil_r+root_i,A,C_l,prec);
                            arb_add(P_each_Y,P_each_Y,w,prec);
                        }
                        
                        arb_set(res,P_each_Y);
                        
                        break;
                    default :
                        printf("PS_Method -> C_l_probability -> Probability_C_l->delta_type->zeta_type 有误\n");
                        exit(1);
                }
                
                _arb_vec_clear(muil_r, root_num);
                arb_clear(Root_Y_parameter->C_l);
                arb_clear(Root_Y_parameter->A);
                free(Root_Y_parameter); //手动释放自定义结构体内存
                
                break;
             default:
                printf("PS_Method -> C_l_probability ->  Probability_C_l -> Power_spectrum_type 有误\n");
                exit(1);
        }
    }
    
    
    
    
    //完成计算，释放
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    arb_clear(A);
    arb_clear(Y);
    arb_clear(t_cl);
    arb_clear(P_each_Y);
    
    return 0;
}


//计算compaction function C 的概率密度分布
int Probability_C(arb_t res, const arb_t C, slong prec)
{
    arb_t s,t,cl;
    
    arb_init(s);
    arb_init(t);
    arb_init(cl);
    
    //P(C)=P(C_l)/(1-3/4*C_l)
    
    //通过C反解C_l
    Trans_C_to_C_l(cl,C,prec);
    
    Probability_C_l(s,cl,prec); //P(C_l)
    
    arb_one(t);
    arb_mul_ui(t,t,3,prec); //1-3/4*C_l
    arb_div_ui(t,t,4,prec);
    arb_mul(t,t,cl,prec);
    arb_neg(t,t);
    arb_add_ui(t,t,1,prec);
    
    arb_div(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(cl);
    
    return 0;
}
