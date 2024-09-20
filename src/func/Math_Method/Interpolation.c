#include "Interpolation.h"
#include <stdlib.h>
#include <arb_poly.h>

//函数插值

//插值点搜索算法: Interpolation search
// https://dl.acm.org/doi/abs/10.1145/359545.359557
//对点数N的平均复杂度为 O(log(log(N)))
//一般的从左到右的遍历 O(N)；二分法（Binary Search） O(log(N))

//关于对 Interpolation search 进一步的优化
//https://lemire.me/blog/2020/11/25/how-fast-does-interpolation-search-converge/
//https://pages.cs.wisc.edu/~chronis/files/efficiently_searching_sorted_arrays.pdf


//位置搜索 
//Interpolation_position_i_search(i,x,x[],N,prec)
slong Interpolation_position_i_search(const arb_t x, const arb_ptr x_array, const slong N, slong prec )
{
    //这里，x[] 的顺序是排好的了，从大到小或从小到大，且值与指标间的关系为均匀分布
    //可以使用 Interpolation search 算法来加速查找
    //这里，能正确查找的前提是数据已排序，对于没有排序的数据，查找会失败，返回 -1
    //查找的数组越均匀分布越快，最慢是O(n)，对纯均匀分布，只需迭代一次即可找到
    //总共 N 个数，指标从 0 开始
    
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    fmpz_t n;
    fmpz_init(n);
    
    slong i;//最终得到的i值，表示 x[i]≤x≤x[i+1]
    slong a,b; //[a,b]
    a=0;
    b=N-1;
    
    int j;
    j=0;
    
    if( arb_lt(x_array+0,x_array+1) ) //升序和降序的处理不同
    {
        //升序排列
        if( arb_lt(x,x_array+0) || arb_gt(x,x_array+N-1) )
        {
            arb_clear(s);
            arb_clear(t);
            fmpz_clear(n);
            return -1;
        }
        
        //对擦边情况进行处理
        if( arb_le(x_array+0,x) && arb_le(x,x_array+1) )
        {
            arb_clear(s);
            arb_clear(t);
            fmpz_clear(n);
            return 0;
        }
        
        if( arb_le(x,x_array+N-1) && arb_le(x_array+N-2,x) )
        {
            arb_clear(s);
            arb_clear(t);
            fmpz_clear(n);
            return N-2;
        }
        
        while(a < b)
        {
            //i=a + trunc[ (b-a)*(x-x[a])/(x[b]-x[a]) ]
            arb_sub(t,x,x_array+a,prec);
            arb_sub(s,x_array+b,x_array+a,prec);
            arb_div(t,t,s,prec);
            arb_mul_si(t,t,b-a,prec);
            arb_trunc(t,t,prec);
            
            arb_get_unique_fmpz(n,t); //得到截断后的整数
            i=fmpz_get_si(n); //上面得到的整数是 fmpz_t 类型，还得再转一下
            i=i+a;
            
            if( arb_lt(x_array+i,x) && arb_lt(x_array+i+1,x) )
            {
                a=i+1;
            }else if ( arb_gt(x_array+i,x) && arb_gt(x_array+i-1,x)  )
            {
                b=i-1;
            }else
            {
                if( arb_le(x_array+i,x) && arb_le(x,x_array+i+1) )
                {
                    ;
                }else if ( arb_ge(x_array+i,x) && arb_ge(x,x_array+i-1) )
                {
                    i=i-1;
                }else
                {
                    //这里可能刚好撞上端点，但前面比较不出来
                    i=i; 
                }
                j=1;
                break;
            }
        }
    }else
    {
        //降序排列
        if( arb_gt(x,x_array+0) || arb_lt(x,x_array+N-1) )
        {
            arb_clear(s);
            arb_clear(t);
            fmpz_clear(n);
            return -1;
        }
        
        //对擦边情况进行处理
        if( arb_ge(x_array+0,x) && arb_ge(x,x_array+1) )
        {
            arb_clear(s);
            arb_clear(t);
            fmpz_clear(n);
            return 0;
        }
        
        if( arb_ge(x,x_array+N-1) && arb_ge(x_array+N-2,x) )
        {
            arb_clear(s);
            arb_clear(t);
            fmpz_clear(n);
            return N-2;
        }
        
        while(a < b)
        {
            //i=a + trunc[ (b-a)*(x-x[a])/(x[b]-x[a]) ]
            arb_sub(t,x,x_array+a,prec);
            arb_sub(s,x_array+b,x_array+a,prec);
            arb_div(t,t,s,prec);
            arb_mul_si(t,t,b-a,prec);
            arb_trunc(t,t,prec);
            
            arb_get_unique_fmpz(n,t); //得到截断后的整数
            i=fmpz_get_si(n); //上面得到的整数是 fmpz_t 类型，还得再转一下
            i=i+a;
            
            if( arb_gt(x_array+i,x) && arb_gt(x_array+i+1,x) )
            {
                a=i+1;
            }else if ( arb_lt(x_array+i,x) && arb_lt(x_array+i-1,x)  )
            {
                b=i-1;
            }else
            {
                if( arb_ge(x_array+i,x) && arb_ge(x,x_array+i+1) )
                {
                    ;
                }else if ( arb_le(x_array+i,x) && arb_le(x,x_array+i-1) )
                {
                    i=i-1;
                }else
                {
                    //这里可能刚好撞上端点，但前面比较不出来
                    i=i; 
                }
                j=1;
                break;
            }
        }
    }
    
    if(j==0) //整个 while 跑完，都没有找到
    {
        i=-1;
    }
    
    arb_clear(s);
    arb_clear(t);
    fmpz_clear(n);
    
    return i;
}



//对于一值插值函数 interp(res,x,x[],y[],coef[],n,prec)
//需要在x[]里找到x的位置 x[i]≤x≤x[i+1],再进行插值运算，运算后的系数保存到 coef[]，下次在相同区间直接调用
//插值的点，保证 x 左右各有三个原函数的值被利用

//这里 coef[] 是一个多项式 arb_poly_t 的数组，需要定义一个函数来进行内存管理
//用来存储多项式系数拟合后得到的结果，相同区间，不用每次计算，加快速度

//初始化插值多项式系数，为其分配内存
// Interp_coe_t 初始化的数目与所求的点数需一致
Interp_coe_t Interpolation_coe_init(slong num)
{
    //通过 Interp_coe_t 声明变量后，需将该结构体指针指向结构体
    Interp_coe_t coe= (Interp_coe_t)calloc(1,sizeof(struct Interpolation_coe_structure));
    
    coe->coe_poly=(arb_poly_t*)calloc(num,sizeof(arb_poly_t));
    coe->num=num;
    
    for(slong i=0; i<num; i++) //对每个系数的多项式表达进行初始化
    {
        arb_poly_init(coe->coe_poly[i]);
    }
    
    return coe;
}


//清理初始化插值多项式系数，释放内存
void Interpolation_coe_clear(Interp_coe_t coe, slong num)
{
    //首先释放每个多项式系数
    for(slong i=0; i<num; i++)
    {
        arb_poly_clear(coe->coe_poly[i]);
    }
    
    free(coe);
}


//通过各个点的函数取值，给出相应的拟合函数 y=func(x)
void Interpolation_fit_func(arb_t res, const arb_t x,
                            const arb_ptr x_array, const arb_ptr y_array, const Interp_coe_t coe, const slong N,
                            slong prec)
{
    //找到 x 在插值点中所处位置
    slong i;
    i=Interpolation_position_i_search(x,x_array,N,prec);
    
    if(i==-1)
    {
        printf("\n所求 x 并不在插值点范围，请检查, N=%li \nx = ", N);
        arb_printn(x, 20,0);printf("\nx[0] = ");
        arb_printn(x_array+0, 20,0);printf("\tx[N-1] = ");
        arb_printn(x_array+N-1, 20,0);printf("\n\n");
        exit(0);
    }
    
    //判断对应系数是否已求解
    if( !arb_poly_is_zero(coe->coe_poly[i]) )
    {
        //已求解，可直接使得前面的结果
        arb_poly_evaluate(res,coe->coe_poly[i],x,prec);
        return;
    }
    
    //未求解
    arb_ptr x_coe,y_coe;
    x_coe=_arb_vec_init(6); //利用 x 前后的六个点来求解
    y_coe=_arb_vec_init(6);
    
    if( i < 3 ) //0，1，2 相同
    {
        arb_set(x_coe+0,x_array+0);
        arb_set(x_coe+1,x_array+1);
        arb_set(x_coe+2,x_array+2);
        arb_set(x_coe+3,x_array+3);
        arb_set(x_coe+4,x_array+4);
        arb_set(x_coe+5,x_array+5);
        
        arb_set(y_coe+0,y_array+0);
        arb_set(y_coe+1,y_array+1);
        arb_set(y_coe+2,y_array+2);
        arb_set(y_coe+3,y_array+3);
        arb_set(y_coe+4,y_array+4);
        arb_set(y_coe+5,y_array+5);
        
    }else if( i > (N-5) ) //  N-4, N-3, N-2, N-1 相同
    {
        arb_set(x_coe+0,x_array+N-6);
        arb_set(x_coe+1,x_array+N-5);
        arb_set(x_coe+2,x_array+N-4);
        arb_set(x_coe+3,x_array+N-3);
        arb_set(x_coe+4,x_array+N-2);
        arb_set(x_coe+5,x_array+N-1);
        
        arb_set(y_coe+0,y_array+N-6);
        arb_set(y_coe+1,y_array+N-5);
        arb_set(y_coe+2,y_array+N-4);
        arb_set(y_coe+3,y_array+N-3);
        arb_set(y_coe+4,y_array+N-2);
        arb_set(y_coe+5,y_array+N-1);
        
    }else
    {
        arb_set(x_coe+0,x_array+i-2);
        arb_set(x_coe+1,x_array+i-1);
        arb_set(x_coe+2,x_array+i);
        arb_set(x_coe+3,x_array+i+1);
        arb_set(x_coe+4,x_array+i+2);
        arb_set(x_coe+5,x_array+i+3);
        
        arb_set(y_coe+0,y_array+i-2);
        arb_set(y_coe+1,y_array+i-1);
        arb_set(y_coe+2,y_array+i);
        arb_set(y_coe+3,y_array+i+1);
        arb_set(y_coe+4,y_array+i+2);
        arb_set(y_coe+5,y_array+i+3);
    }
    
    //利用牛顿插值得到相应的多项式
    
    /***
    //没有复现
    arb_poly_t p_s; //注意，这里用一个临时 arb_poly_t 求解，然后再复制，不要直接使用 coe->coe_poly[i]
                    //因为，在计算中会导致计算的误差范围特别大，使得结果不可用，可能是 arb 的一个 BUG ??
    arb_poly_init(p_s);
    arb_poly_clear(coe->coe_poly[i]);
    arb_poly_init(coe->coe_poly[i]); //还要再重新初始化一下，BUG ??
    
    arb_poly_interpolate_newton(coe->coe_poly[i],x_coe,y_coe,6,prec); //获得插值多项式
    arb_poly_set(coe->coe_poly[i],p_s);
    
    arb_poly_clear(p_s);
    ***/
    
    arb_poly_interpolate_newton(coe->coe_poly[i],x_coe,y_coe,6,prec); //获得插值多项式，牛顿插值法
    arb_poly_evaluate(res,coe->coe_poly[i],x,prec); //利用多项式求解
}


//常微分方程求解中，RFK45 方法对应的内部插值法
void Interpolation_fit_func_odes_RFK45(arb_t res, const arb_t x,
                                       ODEs_RFK45_dense_t dense_out, const slong i_y, //i_y 表示微分方程解 y 中的第 i 个
                                       slong prec)
{
    
    //找到 x 在插值点中所处位置
    slong i,N;
    
    N=dense_out->num_real; //点数个数
    i=Interpolation_position_i_search(x,dense_out->x,N,prec);
    
    if(i==-1)
    {
        printf("\n所求 x 并不在插值点范围，请检查, N=%li \nx = ", N);
        arb_printn(x, 20,0);printf("\nx[0] = ");
        arb_printn(dense_out->x+0, 20,0);printf("\tx[N-1] = ");
        arb_printn(dense_out->x+N-1, 20,0);printf("\n\n");
        exit(0);
    }
    
    arb_t s,t,theta;
    arb_init(s);
    arb_init(t);
    arb_init(theta);
    
    //采用 Hermite interpolation，见 II.6 Dense Output, Discontinuities, Derivatives P190 (6.7)
    //y(xn+θh)=(1-θ)*y_i+θ*y_{i+1}+θ*(θ-1)*( (1-2*θ)(y_{i+1}-y_i)+(θ-1)*h*yp_i+θ*h*yp_{i+1} )
    
    // θ=(x-x_i)/h
    arb_sub(theta,x,dense_out->x+i,prec);
    arb_div(theta,theta,dense_out->h+i,prec);
    
    //y(xn+θh)=(1-θ)*y_0+θ*y_1+θ*(θ-1)*( (1-2*θ)(y_1-y_0)+(θ-1)*h_yp_0+θ*h_yp_1 )
    
    arb_mul_si(t,theta,-2,prec); //(1-2*θ)(y_1-y_0)+(θ-1)*h_yp_0+θ*h_yp_1
    arb_add_si(t,t,1,prec);
    arb_sub(s,dense_out->y_1[i]+i_y,dense_out->y_0[i]+i_y,prec);
    arb_mul(s,s,t,prec);
    
    arb_add_si(t,theta,-1,prec); 
    arb_mul(t,t,dense_out->h_yp_0[i]+i_y,prec);
    arb_add(s,s,t,prec);
    
    arb_mul(t,theta,dense_out->h_yp_1[i]+i_y,prec);
    arb_add(s,s,t,prec);
    
    arb_add_si(t,theta,-1,prec); //θ*(θ-1)
    arb_mul(t,t,theta,prec);
    arb_mul(s,s,t,prec);
    
    
    arb_mul(t,theta,dense_out->y_1[i]+i_y,prec);
    arb_add(s,s,t,prec);
    
    arb_neg(t,theta);
    arb_add_ui(t,t,1,prec);
    arb_mul(t,t,dense_out->y_0[i]+i_y,prec);
    
    arb_add(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(theta);
}



//常微分方程求解中，DOPRI54 方法对应的内部插值法
void Interpolation_fit_func_odes_DOPRI54(arb_t res, const arb_t x,
                                 ODEs_DOPRI54_dense_t dense_out, const slong i_y, //i_y 表示微分方程解 y 中的第 i 个
                                 slong prec)
{
    
    //找到 x 在插值点中所处位置
    slong i,N;
    
    N=dense_out->num_real; //点数个数
    i=Interpolation_position_i_search(x,dense_out->x,N,prec);
    
    if(i==-1)
    {
        printf("\n所求 x 并不在插值点范围，请检查, N=%li \nx = ", N);
        arb_printn(x, 20,0);printf("\nx[0] = ");
        arb_printn(dense_out->x+0, 20,0);printf("\tx[N-1] = ");
        arb_printn(dense_out->x+N-1, 20,0);printf("\n\n");
        exit(0);
    }
    
    arb_t s,t,theta;
    arb_init(s);
    arb_init(t);
    arb_init(theta);
    
    //参看 https://math.stackexchange.com/questions/2947231/how-can-i-derive-the-dense-output-of-ode45
    
    // θ=(x-x_i)/h
    arb_sub(theta,x,dense_out->x+i,prec);
    arb_div(theta,theta,dense_out->h+i,prec);
    
    //y(xn+θh)=r1+θ(r2+(1−θ)(r3+θ(r4+(1−θ)r5)))
    arb_neg(t,theta);
    arb_add_si(t,t,1,prec); //t=(1-θ)
    arb_mul(s,t,dense_out->r_5[i]+i_y,prec);
    
    arb_add(s,s,dense_out->r_4[i]+i_y,prec);
    arb_mul(s,s,theta,prec);
    
    arb_add(s,s,dense_out->r_3[i]+i_y,prec);
    arb_mul(s,s,t,prec);
    
    arb_add(s,s,dense_out->r_2[i]+i_y,prec);
    arb_mul(s,s,theta,prec);
    
    arb_add(res,s,dense_out->r_1[i]+i_y,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(theta);
}


//常微分方程求解中，DOP853 方法对应的内部插值法
void Interpolation_fit_func_odes_DOP853(arb_t res, const arb_t x,
                                        ODEs_DOP853_dense_t dense_out, const slong i_y, //i_y 表示微分方程解 y 中的第 i 个
                                        slong prec)
{
    
    //找到 x 在插值点中所处位置
    slong i,N;
    
    N=dense_out->num_real; //点数个数
    i=Interpolation_position_i_search(x,dense_out->x,N,prec);
    
    if(i==-1)
    {
        printf("\n所求 x 并不在插值点范围，请检查, N=%li \nx = ", N);
        arb_printn(x, 20,0);printf("\nx[0] = ");
        arb_printn(dense_out->x+0, 20,0);printf("\tx[N-1] = ");
        arb_printn(dense_out->x+N-1, 20,0);printf("\n\n");
        exit(0);
    }
    
    arb_t s,s1,t;
    arb_init(s);
    arb_init(s1);
    arb_init(t);
    
    //参看 https://github.com/robclewley/pydstool/blob/master/PyDSTool/integrator/dop853.c
    //s = (x - xold) / hout;
    //s1 = 1.0 - s;
    // rcont1[i]+s*( rcont2[i]+s1*( rcont3[i]+s*( rcont4[i]+s1*( rcont5[i]+s*(rcont6[i]+s1*(rcont7[i]+s*rcont8[i])) ) ) ) );
    
    // s=θ=(x-x_i)/h
    arb_sub(s,x,dense_out->x+i,prec);
    arb_div(s,s,dense_out->h+i,prec);
    
    arb_neg(s1,s); //t1=1-θ=1-s
    arb_add_ui(s1,s1,1,prec);
    
    //λ=rcont5[i]+s*(rcont6[i]+s1*(rcont7[i]+s*rcont8[i]))
    arb_mul(t,s,dense_out->rcont8[i]+i_y,prec);
    arb_add(t,t,dense_out->rcont7[i]+i_y,prec);
    arb_mul(t,t,s1,prec);
    arb_add(t,t,dense_out->rcont6[i]+i_y,prec);
    arb_mul(t,t,s,prec);
    arb_add(t,t,dense_out->rcont5[i]+i_y,prec);
    
    //rcont1[i]+s*( rcont2[i]+s1*( rcont3[i]+s*( rcont4[i]+s1*λ ) ) )
    arb_mul(t,t,s1,prec);
    arb_add(t,t,dense_out->rcont4[i]+i_y,prec);
    arb_mul(t,t,s,prec);
    arb_add(t,t,dense_out->rcont3[i]+i_y,prec);
    arb_mul(t,t,s1,prec);
    arb_add(t,t,dense_out->rcont2[i]+i_y,prec);
    arb_mul(t,t,s,prec);
    arb_add(res,t,dense_out->rcont1[i]+i_y,prec);
    
    
    arb_clear(s);
    arb_clear(s1);
    arb_clear(t);
}

