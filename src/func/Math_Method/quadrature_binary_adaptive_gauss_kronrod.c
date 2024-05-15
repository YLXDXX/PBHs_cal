#include "quadrature_binary_adaptive.h"
#include <stdlib.h>

void small2big_order_2D(arb_ptr x_as, arb_ptr x_bs, arb_ptr y_as, arb_ptr y_bs, arb_ptr x_es, arb_ptr y_es, arb_ptr gs, slong n)
{
    
    //每次新增一个长度，但是有两个是相同的，n-1=n
    //每次排序的工作需要将这两个新数插入整个序列中，形而一个区间单隔从小到大的排列
    //这里个数由少到多，只需要一次冒泡即可
    arb_t inf;
    arb_init(inf);
    arb_pos_inf(inf);
    
    slong i;
    i=n-1;
    
    
    //判断gs是否为无穷大
    if( arb_eq(gs+i,inf) )
    {
        for (; i>0; i--)
        {
            if(arb_lt(x_es+i, x_es+i-1))
            {
                arb_swap(x_as + i, x_as + i-1);
                arb_swap(x_bs + i, x_bs + i-1);
                arb_swap(y_as + i, y_as + i-1);
                arb_swap(y_bs + i, y_bs + i-1);
                arb_swap(x_es + i, x_es + i-1);
                arb_swap(y_es + i, y_es + i-1);
                arb_swap(gs + i, gs + i-1);
                
                arb_swap(x_as + i+1, x_as + i);
                arb_swap(x_bs + i+1, x_bs + i);
                arb_swap(y_as + i+1, y_as + i);
                arb_swap(y_bs + i+1, y_bs + i);
                arb_swap(x_es + i+1, x_es + i);
                arb_swap(y_es + i+1, y_es + i);
                arb_swap(gs + i+1, gs + i);
            }else
            {
                break;
            }
            
        } 
    }else
    {
        for (; i>0; i--)
        {
            if(arb_lt(gs+i, gs+i-1))
            {
                arb_swap(x_as + i, x_as + i-1);
                arb_swap(x_bs + i, x_bs + i-1);
                arb_swap(y_as + i, y_as + i-1);
                arb_swap(y_bs + i, y_bs + i-1);
                arb_swap(x_es + i, x_es + i-1);
                arb_swap(y_es + i, y_es + i-1);
                arb_swap(gs + i, gs + i-1);
                
                arb_swap(x_as + i+1, x_as + i);
                arb_swap(x_bs + i+1, x_bs + i);
                arb_swap(y_as + i+1, y_as + i);
                arb_swap(y_bs + i+1, y_bs + i);
                arb_swap(x_es + i+1, x_es + i);
                arb_swap(y_es + i+1, y_es + i);
                arb_swap(gs + i+1, gs + i);
            }else
            {
                break;
            }
            
        }
        
        
    }
    
    arb_clear(inf);
    
}


//二元函数 Gauss–Kronrod积分法，自适应 对y积分用
static int binary_integration_gauss_kronrod_cal_y(arb_t res, const arb_t f_x,
                                                  my_calc_func_binary func, void *param, const slong order,
                                                  const arb_t a, const arb_t b, const arb_t error, arb_t geterro,
                                                  slong prec) //计算积分
{
    
    arb_t s,t,x,aa,bb,sum_1,sum_2;
    
    arb_init(s);
    arb_init(t);
    arb_init(x);
    arb_init(aa);
    arb_init(bb);
    arb_init(sum_1);
    arb_init(sum_2);
    
    //积分区间[a,b]转到[-1,1]
    //I=(b-a)/2 * w_i*func[(b-a)/2*x_i+(b+a)/2]
    arb_sub(aa,b,a,prec);
    arb_div_si(aa,aa,2,prec);
    
    arb_add(bb,b,a,prec);
    arb_div_si(bb,bb,2,prec);
    
    extern arb_ptr INT_GUSS_KRONROD_COFFI;
    
    //65个点的版本
    for (int i=0; i < 65; i++)
    {
        arb_mul(x,aa,INT_GUSS_KRONROD_COFFI+i,prec);
        arb_add(x,x,bb,prec);
        
        //arb_printn(s, 50,0);printf("\n"); //打印变量
        
        func(s,f_x,x,param,order,prec); //函数值还会再用
        arb_mul(t,s,INT_GUSS_KRONROD_COFFI+65+i,prec);
        
        arb_add(sum_1,sum_1,t,prec);
        
        if( i%2==1 )
        {
            arb_mul(t,s,INT_GUSS_KRONROD_COFFI+129+(i+1)/2,prec);
            
            arb_add(sum_2,sum_2,t,prec);
        }
        
        
    }
    
    arb_mul(sum_1,sum_1,aa,prec);
    arb_mul(sum_2,sum_2,aa,prec);
    
    //误差判定
    arb_set(res,sum_1);
    arb_sub(s,sum_1,sum_2,prec);
    arb_abs(s,s);
    
    arb_set(geterro,s); //返回相应的误差大小
    
    
    if( arb_le(s,error) )
    {
        arb_clear(s);
        arb_clear(t);
        arb_clear(x);
        arb_clear(aa);
        arb_clear(bb);
        arb_clear(sum_1);
        arb_clear(sum_2);
        
        return 0;
    }else{
        
        arb_clear(s);
        arb_clear(t);
        arb_clear(x);
        arb_clear(aa);
        arb_clear(bb);
        arb_clear(sum_1);
        arb_clear(sum_2);
        
        return 1;
    }
    
    
}



//二元函数对y进行积分 自适应计算 迭代版本 iterate

//二元函数 Gauss–Kronrod积分法，自适应 对y积分后再对x积分用
static int binary_integration_gauss_kronrod_cal_y_x(arb_t res,
                                                    my_calc_func_binary func, void *param, const slong order,
                                                    const arb_t y_a, const arb_t y_b, const arb_t y_error,
                                                    const arb_t a, const arb_t b, const arb_t error, arb_t geterro,
                                                    slong prec) //计算积分
{
    
    arb_t s,t,x,aa,bb,sum_1,sum_2,y_geterro,y_geterro_sum;
    
    arb_init(s);
    arb_init(t);
    arb_init(x);
    arb_init(aa);
    arb_init(bb);
    arb_init(sum_1);
    arb_init(sum_2);
    arb_init(y_geterro);
    arb_init(y_geterro_sum);
    
    //y计算的误差初始为零
    arb_zero(y_geterro_sum);
    
    //初始求和为零
    arb_zero(sum_1);
    arb_zero(sum_2);
    
    //积分区间[a,b]转到[-1,1]
    //I=(b-a)/2 * w_i*func[(b-a)/2*x_i+(b+a)/2]
    arb_sub(aa,b,a,prec);
    arb_div_si(aa,aa,2,prec);
    
    arb_add(bb,b,a,prec);
    arb_div_si(bb,bb,2,prec);
    
    extern arb_ptr INT_GUSS_KRONROD_COFFI;
    
    
    //65个点的版本
    for (int i=0; i < 65; i++)
    {
        arb_mul(x,aa,INT_GUSS_KRONROD_COFFI+i,prec);
        arb_add(x,x,bb,prec);
        
        //arb_printn(s, 50,0);printf("\n"); //打印变量
        
        //此时，计算的函数值是对y积分后的结果
        //func(s,f_x,x,param,order,prec); //函数值还会再用
        
        binary_integration_gauss_kronrod_cal_y(s,x,func,param,order,y_a,y_b,y_error,y_geterro,prec);
        
        arb_add(y_geterro_sum,y_geterro_sum,y_geterro,prec);
        
        arb_mul(t,s,INT_GUSS_KRONROD_COFFI+65+i,prec);
        
        arb_add(sum_1,sum_1,t,prec);
        
        if( i%2==1 )
        {
            arb_mul(t,s,INT_GUSS_KRONROD_COFFI+129+(i+1)/2,prec);
            
            arb_add(sum_2,sum_2,t,prec);
        }
    }
    
    arb_mul(sum_1,sum_1,aa,prec);
    arb_mul(sum_2,sum_2,aa,prec);
    
    arb_set(res,sum_1);
    
    //x误差判定
    arb_sub(s,sum_1,sum_2,prec);
    arb_abs(s,s);
    
    
    //arb_set(geterro,s); //返回相应的误差大小
    arb_div_ui(y_geterro_sum,y_geterro_sum,65,prec); //y误差求平均
    arb_add(geterro,s,y_geterro_sum,prec); //x误差加y误差
    
    int res_lable;
    if( arb_le(s,error) && arb_le(y_geterro_sum,y_error) ) //x和y均满足
    {
        res_lable=00;
    }else if ( arb_le(s,error) ) //x满足
    {
        res_lable=01;
    } else if ( arb_le(y_geterro_sum,y_error) ) //y满足
    {
        res_lable=10;
    } else //x和y均不满足
    {
        res_lable=11;
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(x);
    arb_clear(aa);
    arb_clear(bb);
    arb_clear(sum_1);
    arb_clear(sum_2);
    arb_clear(y_geterro);
    arb_clear(y_geterro_sum);
    
    return res_lable;
}



//自适应计算 迭代版本 iterate
int integration_binary_rectangle_adaptive_gauss_kronrod(arb_t res, my_calc_func_binary func, void *param, const slong order,
                                          const arb_t x_a, const arb_t x_b, const arb_t x_error,
                                          slong x_step_min , slong x_step_max,
                                          const arb_t y_a, const arb_t y_b, const arb_t y_error,
                                          slong y_step_min , slong y_step_max, 
                                          slong prec)
{
    get_gauss_kronrod_node_weight(65,prec); //获取 gauss_kronrod 的节点位置和权重
    arb_ptr x_as, x_bs, y_as, y_bs, x_es, y_es, gs;
    arb_t s, u, w;
    
    slong depth,top,alloc,leaf_interval_count;
    int judge,stopping, status;
    slong step_max,step_min;
    
    step_min=x_step_min*y_step_min;
    step_max=x_step_max*y_step_max;
    
    judge=11; //默认x和y均分割
    depth=1;
    stopping = 0;
    status=0;
    alloc = 16;
    leaf_interval_count=0;
    
    x_as = _arb_vec_init(alloc); //初始化
    x_bs = _arb_vec_init(alloc);
    y_as = _arb_vec_init(alloc);
    y_bs = _arb_vec_init(alloc);
    x_es = _arb_vec_init(alloc); //参考误差
    y_es = _arb_vec_init(alloc);
    gs = _arb_vec_init(alloc); //计算得到的实际误差
    
    arb_set(x_as, x_a);//初始积分区间
    arb_set(x_bs, x_b);
    arb_set(y_as, y_a);
    arb_set(y_bs, y_b);
    arb_set(x_es, x_error); //用来存放计算所需对比的误差
    arb_set(y_es, y_error);
    arb_pos_inf(gs); //初始化计算误差
    
    arb_init(s);//积分求和用
    arb_init(u);
    arb_init(w);
    
    arb_pos_inf(w);//最初设其为正无穷大，用于最小迭代次数内设置计算误差
    arb_zero(u); //最初设为零，用于最小迭代次数内设置计算结果
    arb_zero(s);
    
    int time_min;
    int min_interval;
    time_min=0; //判断最小迭代次数
    min_interval=0; //判断是否达到最小区间间隔
    
    int ret_judge=0;
    
    arb_t temp_t; //临时变量
    arb_init(temp_t);
    
    while (depth >= 1)
    {
        
        //超过最大迭代次数，终止
        if (stopping == 0 && depth >= step_max )
        {
            //flint_printf("\n达到最大迭代深度：  %wd\n\n", step_max);
            status++;
            stopping = 1;
            
            ret_judge=1;
            
            /*
             *           for(slong i=0;i<=depth;i++)
             *           {
             *               //arb_sub(u,bs+i,as+i,prec);
             *               arb_printn(gs+i, 15,0);printf("\n");
        }
        exit(0);
        */
        }
        
        
        top = depth - 1;
        
        if ( stopping) //达到最大迭代后，结束
        {
            //直接计算每个区间的值，加起来
            binary_integration_gauss_kronrod_cal_y_x(u,func,param,order,
                                                     y_as+top,y_bs+top,y_es+top,
                                                     x_as+top,x_bs+top,x_es+top,w,prec);
            
            arb_add(s, s, u, prec);
            
            leaf_interval_count++;
            
            depth--;
            
            continue;
        }
        
        //有最少的初始分隔次数要求
        if( time_min==0  && depth >= step_min )
        {
            time_min=1;
        }
        
        if( time_min ) //有最少积分区间分隔次数要求
        {
            judge = binary_integration_gauss_kronrod_cal_y_x(u,func,param,order,
                                                             y_as+top,y_bs+top,y_es+top,
                                                             x_as+top,x_bs+top,x_es+top,w,prec);
            
            if( judge==00 || min_interval==1 ) // min_interval==1 用于判断是否二分到最小区间间隔
            {
                //满足精度要求
                arb_add(s, s, u, prec);
                
                leaf_interval_count++;
                
                depth--;
                
                continue;
            }
        }
        
        //看数组存储是否够用，不够再增加分配
        if (depth >= alloc - 3)
        {
            slong k;
            x_as = flint_realloc(x_as, 2 * alloc * sizeof(arb_struct));
            x_bs = flint_realloc(x_bs, 2 * alloc * sizeof(arb_struct));
            y_as = flint_realloc(y_as, 2 * alloc * sizeof(arb_struct));
            y_bs = flint_realloc(y_bs, 2 * alloc * sizeof(arb_struct));
            x_es = flint_realloc(x_es, 2 * alloc * sizeof(arb_struct));
            y_es = flint_realloc(y_es, 2 * alloc * sizeof(arb_struct));
            gs = flint_realloc(gs, 2 * alloc * sizeof(arb_struct));
            for (k = alloc; k < 2 * alloc; k++)
            {
                arb_init(x_as + k); //分配空间后初始化
                arb_init(x_bs + k);
                arb_init(y_as + k);
                arb_init(y_bs + k);
                arb_init(x_es + k);
                arb_init(y_es + k);
                arb_init(gs + k);
            }
            alloc *= 2;
        }
        
        
        //根据 judge=00/11/10/01 判断分割区间
        if ( judge==11 ) //x和y都需要分割
        {
            //给定一个矩形，按中点分成四个矩形
            //二分区间 [x_a,x_mid],[x_mid,x_b],[y_a,y_mid],[y_mid,y_b]
            //有四个积分区域：[x_a,x_mid,y_a,y_mid], [x_a,x_mid,y_mid,y_b],[x_mid,x_b,y_a,y_mid],[x_mid,x_b,y_mid,y_b]
            //x_as --> x_a,  x_a,  x_mid,x_mid
            //x_bs --> x_mid,x_mid,x_b,  x_b
            //y_as --> y_a,  y_mid,y_a,  y_mid
            //y_bs --> y_mid,y_b,  y_mid,y_b
            
            arb_set(x_as + top+1, x_as+top);
            arb_set(x_bs + top+2, x_bs+top);
            arb_set(x_bs + top+3, x_bs+top);
            
            arb_set(y_as + top+2, y_as+top);
            arb_set(y_bs + top+1, y_bs+top);
            arb_set(y_bs + top+3, y_bs+top);
            
            arb_add(x_as + top+2, x_as+top, x_bs+top, prec);
            arb_mul_2exp_si(x_as + top+2, x_as + top+2, -1);
            arb_set(x_as + top+3,x_as + top+2);
            arb_set(x_bs + top,x_as + top+2);
            arb_set(x_bs + top+1,x_as + top+2);
            
            arb_add(y_as + top+1, y_as+top, y_bs+top, prec);
            arb_mul_2exp_si(y_as + top+1, y_as + top+1, -1);
            arb_set(y_as + top+3,y_as + top+1);
            arb_set(y_bs + top,y_as + top+1);
            arb_set(y_bs + top+2,y_as + top+1);
            
            
            //区间误差更新
            arb_mul_2exp_si(x_es+top, x_es+top, -2); // es/4
            arb_set(x_es + depth , x_es+top); // es/4
            arb_set(x_es + top+2 , x_es+top);
            arb_set(x_es + top+3 , x_es+top);
            
            
            arb_mul_2exp_si(y_es+top, y_es+top, -2); // es/4
            arb_set(y_es + depth , y_es+top); // es/4
            arb_set(y_es + top+2 , y_es+top);
            arb_set(y_es + top+3 , y_es+top);
            
            //计算得到的区间误差更新
            arb_mul_2exp_si(gs+depth, w, -2); //暴力处理，将原区间计算误差除以4        
            arb_set(gs+top, gs+depth);
            arb_set(gs+depth+1, gs+depth);
            arb_set(gs+depth+2, gs+depth);
            
            //判断是否达到最小子区间间隔,这里只判断了x的区间
            arb_sub(temp_t,x_bs + depth+2,x_as + depth+2,prec);
            arb_abs(temp_t,temp_t); //求绝对值，积分区间大小可能相反
            
            depth=depth+3; //每种情况不同
            
        } else if ( judge==10 ) //x需要分割，y不需要
        {
            //两个积分区域 [x_a,x_mid,y_a,y_b] [x_mid,x_b,y_a,y_b]
            //x_as --> x_a,  x_mid
            //x_bs --> x_mid,x_b
            //y_as --> y_a,  y_a
            //y_bs --> y_b,  y_b
            
            arb_set(x_bs + top+1, x_bs+top);
            arb_set(y_as + top+1, y_as+top);
            arb_set(y_bs + top+1, y_bs+top);
            
            arb_add(x_as + top+1, x_as+top, x_bs+top, prec);
            arb_mul_2exp_si(x_as + top+1, x_as + top+1, -1);
            arb_set(x_bs+top, x_as + top+1);
            
            //区间误差更新
            arb_mul_2exp_si(x_es+top, x_es+top, -1); // es/2
            arb_set(x_es + depth , x_es+top); // es/2
            
            arb_mul_2exp_si(y_es+top, y_es+top, -1); // es/2
            arb_set(y_es + depth , y_es+top); // es/2
            
            //计算得到的区间误差更新
            arb_mul_2exp_si(gs+depth, w, -1); //暴力处理，将原区间计算误差除以2        
            arb_set(gs+top, gs+depth);
            
            //判断是否达到最小子区间间隔,这里只判断了x的区间
            arb_sub(temp_t,x_bs + depth,x_as + depth,prec);
            arb_abs(temp_t,temp_t); //求绝对值，积分区间大小可能相反
            
            depth=depth+1; //每种情况不同
            
        }else //judge==01 //y需要分割，x不需要
        {
            //两个积分区域 [x_a,x_b,y_a,y_mid] [x_a,x_b,y_mid,y_b]
            //x_as --> x_a,  x_a
            //x_bs --> x_b,  x_b
            //y_as --> y_a,  y_mid
            //y_bs --> y_mid,y_b
            
            arb_set(x_as + top+1, x_as+top);
            arb_set(x_bs + top+1, x_bs+top);
            arb_set(y_bs + top+1, y_bs+top);
            
            arb_add(y_as + top+1, y_as+top, y_bs+top, prec);
            arb_mul_2exp_si(y_as + top+1, y_as + top+1, -1);
            arb_set(y_bs+top, y_as + top+1);
            
            //区间误差更新
            arb_mul_2exp_si(x_es+top, x_es+top, -1); // es/2
            arb_set(x_es + depth , x_es+top); // es/2
            
            arb_mul_2exp_si(y_es+top, y_es+top, -1); // es/2
            arb_set(y_es + depth , y_es+top); // es/2
            
            //计算得到的区间误差更新
            arb_mul_2exp_si(gs+depth, w, -1); //暴力处理，将原区间计算误差除以2        
            arb_set(gs+top, gs+depth);
            
            //判断是否达到最小子区间间隔,这里只判断了x的区间
            arb_sub(temp_t,y_bs + depth,y_as + depth,prec);
            arb_abs(temp_t,temp_t); //求绝对值，积分区间大小可能相反
            
            depth=depth+1; //每种情况不同
            
        }
        
        //判断是否达到最小子区间间隔
        if(arb_lt(temp_t,INT_MIN_INTERVAL)) //对于一些间断点，区间可能会无限分隔下去，需设定最小区间间隔
        {
            min_interval=1; //达到最小区间间隔，强制积分求和
        }else
        {
            min_interval=0;
        }
        
        
        if(1)
        {
            //区间单隔，大的上升，小的沉底，按从小到大的顺序排列
            small2big_order_2D(x_as, x_bs, y_as, y_bs, x_es, y_es, gs, depth);
        }
        
        //depth=depth+3;
        //printf("deptht: %li\n",depth);
    }
    
    arb_set(res, s);
    
    //printf("leaf_interval_count: %li\n",leaf_interval_count);
    
    _arb_vec_clear(x_as, alloc);//释放内存
    _arb_vec_clear(x_bs, alloc);
    _arb_vec_clear(y_as, alloc);
    _arb_vec_clear(y_bs, alloc);
    _arb_vec_clear(x_es, alloc);
    _arb_vec_clear(y_es, alloc);
    _arb_vec_clear(gs, alloc);
    arb_clear(s);
    arb_clear(u);
    arb_clear(w);
    arb_clear(temp_t);
    
    if(ret_judge==0)
    {
        return 0;
    }else
    {
        return 1;
    }
}


