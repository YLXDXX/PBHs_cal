#include "quadrature.h" 
#include <stdlib.h>
#include <omp.h>

//自适应辛普森积分
//可以指定精度
//Simpson公式
int simpson_formula(arb_t res, int (*func)(arb_t f_res, const arb_t x, slong prec),
                    const arb_t a, const arb_t b, slong prec)
{
    //局部变量
    arb_t mid,t;
    
    arb_init(mid);
    arb_init(t);
    slong n;
    
    //(f(a) + 4.0 * f(mid) + f(b) ) * (b - a) / 6.0;
    
    arb_add(mid,a,b,prec);
    n=2;
    arb_div_si(mid,mid,n,prec);
    
    func(res,a,prec);
    
    func(t,mid,prec);
    n=4;
    arb_mul_si(t,t,n,prec);
    arb_add(res,res,t,prec);
    
    func(t,b,prec);
    arb_add(res,res,t,prec);
    
    arb_sub(t,b,a,prec);
    arb_mul(res,res,t,prec);
    
    n=6;
    arb_div_si(res,res,n,prec);
    
    //完成计算，释放
    arb_clear(mid);
    arb_clear(t);
    return 0;
}

//Simpsont积分，加入最少迭代次数，加入最大迭代次数
int integration_simpson(arb_t int_sum, int (*func)(arb_t f_res, const arb_t x, slong prec),
                        const arb_t a, const arb_t b,const arb_t ans,const arb_t eps,
                        slong step_min , slong step_max, slong prec)
{
    //局部变量
    arb_t mid,t,s,ll,rr,e;
    
    arb_init(mid);
    arb_init(t);
    arb_init(s);
    arb_init(ll);
    arb_init(rr);
    arb_init(e);
    slong n;
    
    
    arb_add(mid,a,b,prec);
    n=2;
    arb_div_si(mid,mid,n,prec);
    
    //arb_printn(a, 10,0); flint_printf("\n");
    //arb_printn(b, 10,0); flint_printf("\n");
    //arb_printn(mid, 10,0); flint_printf("\n");flint_printf("\n");flint_printf("\n");
    
    //ll = simpson(a, mid), rr = simpson(mid, b);
    
    simpson_formula(ll,func,a,mid,prec);
    simpson_formula(rr,func,mid,b,prec);
    
    //确认精度
    //if (fabs(ll + rr - ans) <= 15 * eps) 
    arb_add(t,ll,rr,prec);
    arb_sub(t,t,ans,prec);
    arb_abs(t,t); //绝对值
    n=15;
    arb_mul_si(s,eps,n,prec);
    
    arb_sub(t,s,t,prec);
    
    //arb_printn(mid, 30, 0);flint_printf("22\n");
    
    
    
    if ( arb_is_positive(t) || step_max <= 0 ) 
    {
        
        //有最少迭代次数要求
        if (step_min <= 0)
        {
            //printf("%ld\n",step_min);
            //达到精度要求
            //ll + rr + (ll + rr - ans) / 15.0;
            arb_add(t,ll,rr,prec);
            arb_sub(int_sum,t,ans,prec);
            n=15;
            arb_div_si(int_sum,int_sum,n,prec);
            arb_add(int_sum,int_sum,t,prec);
        }else
        {
            //递归调用，提高精度
            n=2;
            arb_div_si(e,eps,n,prec);
            
            integration_simpson(t, func,a,mid,ll,e,step_min-1,step_max-1,prec);
            integration_simpson(s, func,mid,b,rr,e,step_min-1,step_max-1,prec);
            arb_add(int_sum,t,s,prec);
        }
        
    }else if ( arb_is_nonpositive(t) ) 
    {
        //递归调用，提高精度
        n=2;
        arb_div_si(e,eps,n,prec);
        
        integration_simpson(t, func,a,mid,ll,e,step_min-1,step_max-1,prec);
        integration_simpson(s, func,mid,b,rr,e,step_min-1,step_max-1,prec);
        arb_add(int_sum,t,s,prec);
        
    }else{
        ;
    }
    
    //完成计算，释放
    arb_clear(mid);
    arb_clear(t);
    arb_clear(s);
    arb_clear(ll);
    arb_clear(rr);
    arb_clear(e);
    return 0;
}



//Gauss–Kronrod积分法，自适应，15点的
int integration_gauss_kronrod_cal(arb_t res, my_calc_func func, void *param, const slong order,
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
    
    //初始求和为零
    arb_zero(sum_1);
    arb_zero(sum_2);
    
    //积分区间[a,b]转到[-1,1]
    //I=(b-a)/2 * w_i*func[(b-a)/2*x_i+(b+a)/2]
    arb_sub(aa,b,a,prec);
    arb_div_si(aa,aa,2,prec);
    
    arb_add(bb,b,a,prec);
    arb_div_si(bb,bb,2,prec);
    
    //extern arb_ptr INT_GUSS_KRONROD_COFFI;
    
    //15个点的版本
    /*
    for (int i=0; i < 15; i++)
    {
        arb_mul(x,aa,INT_GUSS_KRONROD_COFFI+i,prec);
        arb_add(x,x,bb,prec);
        
        //arb_printn(s, 50,0);printf("\n"); //打印变量
        
        func(s,x,param,order,prec); //函数值还会再用
        arb_mul(t,s,INT_GUSS_KRONROD_COFFI+15+i,prec);
        arb_mul(t,t,aa,prec);
        arb_add(sum_1,sum_1,t,prec);
        
        if( i%2==1 )
        {
            arb_mul(t,s,INT_GUSS_KRONROD_COFFI+29+(i+1)/2,prec);
            arb_mul(t,t,aa,prec);
            arb_add(sum_2,sum_2,t,prec);
        }
        
        
    }
    */
    
    //65个点的版本
    for (int i=0; i < 65; i++)
    {
        arb_mul(x,aa,INT_GUSS_KRONROD_COFFI+i,prec);
        arb_add(x,x,bb,prec);
        
        //arb_printn(s, 50,0);printf("\n"); //打印变量
        
        func(s,x,param,order,prec); //函数值还会再用
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


void small2big_order(arb_ptr as, arb_ptr bs, arb_ptr es, arb_ptr gs, slong n)
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
            if(arb_lt(es+i, es+i-1))
            {
                arb_swap(as + i, as + i-1);
                arb_swap(bs + i, bs + i-1);
                arb_swap(es + i, es + i-1);
                arb_swap(gs + i, gs + i-1);
                
                arb_swap(as + i+1, as + i);
                arb_swap(bs + i+1, bs + i);
                arb_swap(es + i+1, es + i);
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
                arb_swap(as + i, as + i-1);
                arb_swap(bs + i, bs + i-1);
                arb_swap(es + i, es + i-1);
                arb_swap(gs + i, gs + i-1);
                
                arb_swap(as + i+1, as + i);
                arb_swap(bs + i+1, bs + i);
                arb_swap(es + i+1, es + i);
                arb_swap(gs + i+1, gs + i);
            }else
            {
                break;
            }
            
        }
        
        
    }
    
    arb_clear(inf);
    
}




//自适应计算 迭代版本 iterate
int integration_gauss_kronrod_iterate(arb_t res, my_calc_func func, void *param, const slong order,
                              const arb_t a, const arb_t b, const arb_t error,
                              slong step_min , slong step_max,
                              slong prec)
{
    arb_ptr as, bs, es,gs;
    arb_t s, u, w;
    
    slong depth,top,alloc,leaf_interval_count;
    int judge,stopping, status;
    
    
    depth=1;
    stopping = 0;
    status=0;
    alloc = 4;
    leaf_interval_count=0;
    
    as = _arb_vec_init(alloc); //初始化
    bs = _arb_vec_init(alloc);
    es = _arb_vec_init(alloc); //参考误差
    gs = _arb_vec_init(alloc); //计算得到的实际误差
    
    arb_set(as, a);//初始积分区间
    arb_set(bs, b);
    arb_set(es, error); //用来存放计算所需对比的误差
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
            for(slong i=0;i<=depth;i++)
            {
                //arb_sub(u,bs+i,as+i,prec);
                arb_printn(gs+i, 15,0);printf("\n");
            }
            exit(0);
            */
        }
        
        
        top = depth - 1;
        
        if ( stopping) //达到最大迭代后，结束
        {
            //直接计算每个区间的值，加起来
            integration_gauss_kronrod_cal(u,func,param,order,as+top,bs+top,es+top,w,prec);
            
            arb_add(s, s, u, prec);
            
            leaf_interval_count++;
            
            depth--;
            
            continue;
        }
        
        //有最少的初始分隔次数要求
        if( time_min==0  && depth >= step_min   )
        {
            time_min=1;
        }
            
        if( time_min ) //有最少积分区间分隔次数要求
        {
            judge = integration_gauss_kronrod_cal(u,func,param,order,as+top,bs+top,es+top,w,prec);
            
            if( judge==0 || min_interval==1 ) // min_interval==1 用于判断是否二分到最小区间间隔
            {
                //满足精度要求
                arb_add(s, s, u, prec);
                
                leaf_interval_count++;
                
                depth--;
                
                continue;
            }
        }
        
        //看数组存储是否够用，不够再增加分配
        if (depth >= alloc - 1)
        {
            slong k;
            as = flint_realloc(as, 2 * alloc * sizeof(arb_struct));
            bs = flint_realloc(bs, 2 * alloc * sizeof(arb_struct));
            es = flint_realloc(es, 2 * alloc * sizeof(arb_struct));
            gs = flint_realloc(gs, 2 * alloc * sizeof(arb_struct));
            for (k = alloc; k < 2 * alloc; k++)
            {
                arb_init(as + k); //分配空间后初始化
                arb_init(bs + k);
                arb_init(es + k);
                arb_init(gs + k);
            }
            alloc *= 2;
        }
        
        //二分区间
        /* Interval [depth] becomes [mid, b]. */
        arb_set(bs + depth, bs+top);
        arb_add(as + depth, as+top, bs+top, prec);
        arb_mul_2exp_si(as + depth, as + depth, -1);
        
        /* Interval [0] becomes [a, mid]. */
        arb_set(bs+top, as + depth);
        
        //判断是否达到最小子区间间隔
        arb_sub(temp_t,bs + depth,as + depth,prec);
        arb_abs(temp_t,temp_t); //求绝对值，积分区间大小可能相反
        if(arb_lt(temp_t,INT_MIN_INTERVAL)) //对于一些间断点，区间可能会无限分隔下去，需设定最小区间间隔
        {
           min_interval=1; //达到最小区间间隔，强制积分求和
        }else
        {
            min_interval=0;
        }
        
        //区间误差更新
        arb_mul_2exp_si(es+top, es+top, -1); // es/2
        arb_set(es + depth , es+top); // es/2
        
        
        //计算得到的区间误差更新
        arb_mul_2exp_si(gs+depth, w, -1); //暴力处理，将原区间计算误差除以2        
        arb_set(gs+top, gs+depth);
        
        if(1)
        {
        //区间单隔，大的上升，小的沉底，按从小到大的顺序排列
        small2big_order(as, bs, es, gs, depth);
        }
        //printf("depth: %li\n",depth);
        
        depth++;
        //printf("deptht: %li\n",depth);
    }
    
    arb_set(res, s);
    
    //printf("leaf_interval_count: %li\n",leaf_interval_count);
    
    _arb_vec_clear(as, alloc);//释放内存
    _arb_vec_clear(bs, alloc);
    _arb_vec_clear(es, alloc);
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



//自适应计算 递归版本 recursive
int integration_gauss_kronrod_recursive(arb_t res, my_calc_func func, void *param, const slong order,
                              const arb_t a, const arb_t b, const arb_t error,
                              slong step_min , slong step_max,
                              slong prec)
{
    int judge;
    arb_t aa,bb,s,t,mid,e,w;
    
    arb_init(aa);
    arb_init(bb);
    arb_init(s);
    arb_init(t);
    arb_init(mid);
    arb_init(e);
    arb_init(w);
    
    arb_set(aa,a);
    arb_set(bb,b);
    
    judge = integration_gauss_kronrod_cal(s,func,param,order,aa,bb,error,w, prec);
    
    if(judge==0 || step_max <= 0 ){
        if ( step_min < 0 ){
            //满足精度要求
            arb_set(res,s);
        }else{
            //不满足精度要求
            arb_add(mid,bb,aa,prec);
            arb_div_si(mid,mid,2,prec);
            arb_div_si(e,error,2,prec);//在子区间上，误差也要折半，才以保证总精度
            
            integration_gauss_kronrod_recursive(s,func,param,order,aa,mid,e,step_min-1, step_max-1,prec);
            integration_gauss_kronrod_recursive(t,func,param,order,mid,bb,e,step_min-1, step_max-1,prec);
            
            arb_add(res,s,t,prec);
        }
        
    }else{
        //不满足精度要求
        arb_add(mid,bb,aa,prec);
        arb_div_si(mid,mid,2,prec);
        arb_div_si(e,error,2,prec);//在子区间上，误差也要折半，才以保证总精度
        
        integration_gauss_kronrod_recursive(s,func,param,order,aa,mid,e,step_min-1, step_max-1,prec);
        integration_gauss_kronrod_recursive(t,func,param,order,mid,bb,e,step_min-1, step_max-1,prec);
        
        arb_add(res,s,t,prec);
    }
    
    arb_clear(aa);
    arb_clear(bb);
    arb_clear(s);
    arb_clear(t);
    arb_clear(mid);
    arb_clear(e);
    arb_clear(w);
    
    return 0;
}



// return optimized Exp-Sinh integral split point d
static void exp_sinh_opt_d(arb_t res, my_calc_func func, void *param, const slong order,
                           const arb_t a, const arb_t eps, const arb_t dd, slong prec)
{
    arb_t h2,r,fl,fr,h,s,lfl,lfr,lr,d,t_t,t_s,t_w;
    arb_init(h2);
    arb_init(r);
    arb_init(fl);
    arb_init(fr);
    arb_init(h);
    arb_init(s);
    arb_init(lfl);
    arb_init(lfr);
    arb_init(lr);
    arb_init(d);
    arb_init(t_t);
    arb_init(t_s);
    arb_init(t_w);
    
    arb_set(d,dd); //d=dd
    
    //h2 = f(a + d/2) - f(a + d*2)*4
    arb_mul_ui(t_t,d,2,prec);//f(a + d*2)*4
    arb_add(t_t,t_t,a,prec);
    func(t_w,t_t,param,order,prec); 
    arb_mul_ui(t_w,t_w,4,prec);
    
    arb_mul_2exp_si(t_s,d,-1); //f(a + d/2)
    arb_add(t_s,t_s,a,prec);
    func(t_t,t_s,param,order,prec);
    
    arb_sub(h2,t_t,t_w,prec);
    
    slong i=1,j=32; // j=32 is optimal to find r
    
    arb_abs(t_s,h2);
    arb_one(t_t);
    arb_mul_2exp_si(t_t,t_t,-16); //2^(-16)
    //arb_set_str(t_t,"1e-5",prec);
    if( arb_is_finite(h2) && arb_gt(t_s,t_t) ) //isfinite(h2) && fabs(h2) > 2^-16 // if |h2| > 2^-16
    {
        arb_zero(s); //s=0
        arb_one(lr); //lr=2
        arb_mul_ui(lr,lr,2,prec);
        do{  // find max j such that fl and fr are finite 
            j/=2;
            
            arb_one(r); //r = 1 << (i + j) = 1*2^(i+j)
            arb_mul_2exp_si(r,r,i+j);
            
            arb_div(t_s,d,r,prec); //fl = f(a + d/r)
            arb_add(t_s,t_s,a,prec);
            func(fl,t_s,param,order,prec);
            
            arb_mul(t_s,d,r,prec); //fr = f(a + d*r)*r*r
            arb_add(t_s,t_s,a,prec);
            func(fr,t_s,param,order,prec);
            arb_sqr(t_s,r,prec);
            arb_mul(fr,fr,t_s,prec);
            
            arb_sub(h,fl,fr,prec); // h = fl - fr
            
        } while ( j>1 && !arb_is_finite(h) );
        
        if( j>1 && arb_is_finite(h) && arb_sgn_nonzero(h) != arb_sgn_nonzero(h2) )
        {
            arb_set(lfl,fl); // last fl=f(a+d/r)
            arb_set(lfr,fr); // last fr=f(a+d*r)*r*r 
            do{ // bisect in 4 iterations 
                j/=2;
                
                arb_one(r); //r = 1 << (i + j) = 1*2^(i+j)
                arb_mul_2exp_si(r,r,i+j);
                
                arb_div(t_s,d,r,prec); //fl = f(a + d/r)
                arb_add(t_s,t_s,a,prec);
                func(fl,t_s,param,order,prec);
                
                arb_mul(t_s,d,r,prec); //fr = f(a + d*r)*r*r
                arb_add(t_s,t_s,a,prec);
                func(fr,t_s,param,order,prec);
                arb_sqr(t_s,r,prec);
                arb_mul(fr,fr,t_s,prec);
                
                arb_sub(h,fl,fr,prec); // h = fl - fr
                
                if(arb_is_finite(h))
                {
                    arb_abs(t_s,h); //s += fabs(h) // sum |h| to remove noisy cases
                    arb_add(s,s,t_s,prec);
                    
                    if( arb_sgn_nonzero(h)==arb_sgn_nonzero(h2) )
                    { // search right half
                        i+=j; 
                    }else // search left half
                    {
                        arb_set(lfl,fl); // record last fl=f(a+d/r)
                        arb_set(lfr,fr); // record last fl=f(a+d*r)*r*r 
                        arb_set(lr,r); // record last r
                    }
                    
                }
                
            }while(j>1);
            
            if(arb_gt(s,eps)) //s>eps // if sum of |h| > eps
            {
                arb_sub(h,lfl,lfr,prec); //h = lfl - lfr // use last fl and fr before the sign change 
                arb_set(r,lr);
                if( !arb_is_zero(h) ) //h!=0 // if last diff != 0, back up r by one step
                {
                    arb_mul_2exp_si(r,r,-1); //r/=2
                }
                
                arb_abs(t_s,lfl);
                arb_abs(t_t,lfr);
                if(arb_lt(t_s,t_t))
                {
                    arb_div(d,d,r,prec); // move d closer to the finite endpoint
                }else
                {
                    arb_mul(d,d,r,prec); // move d closer to the infinite endpoint 
                }
                
            }
        }
    }
    
    //arb_set(res,d);
    arb_get_mid_arb(res,d);
    
    arb_clear(h2);
    arb_clear(r);
    arb_clear(fl);
    arb_clear(fr);
    arb_clear(h);
    arb_clear(s);
    arb_clear(lfl);
    arb_clear(lfr);
    arb_clear(lr);
    arb_clear(d);
    arb_clear(t_t);
    arb_clear(t_s);
    arb_clear(t_w);
}


//Double Exponential积分，可适用于无穷区间积分，及区间端点处函数发散的情况
//特别是在区间端点处函数发散的情况，Double Exponential积分明显好于Gauss–Kronrod积分
//From： https://www.genivia.com/qthsh.html
//the improved quad routine with separate branches optimized for Tanh-Sinh, Exp-Sinh and Sinh-Sinh
//下面是其 arb 版本的实现

int Double_Exponential_Quadrature(arb_t res, my_calc_func func, void *param, const slong order,
                                  const arb_t a, const arb_t b, const arb_t eps,
                                  slong step_min , slong step_max,
                                  slong prec)
{
    //函数f(x)在[a,b]上进行积分，区间端点a或b都可以为 ∞，甚至在[-∞，+∞]上积分
    //[a,b] --> [aa,bb]
    arb_t tol,c,d,s,e,v,h,aa,bb,p,q,fp,fm,t,eh,u,r,w,x,y,t_s,t_t,t_w; 
    slong sign=1;
    
    arb_init(tol);
    arb_init(c);
    arb_init(d);
    arb_init(s);
    arb_init(e);
    arb_init(v);
    arb_init(h);
    arb_init(aa);
    arb_init(bb);
    arb_init(p);
    arb_init(q);
    arb_init(fp);
    arb_init(fm);
    arb_init(t);
    arb_init(eh);
    arb_init(u);
    arb_init(r);
    arb_init(w);
    arb_init(x);
    arb_init(y);
    arb_init(t_s);
    arb_init(t_t);
    arb_init(t_w);
    
    //设定初始条件
    arb_set(aa,a);
    arb_set(bb,b);
    arb_mul_ui(tol,eps,10,prec); //tol=10*eps
    arb_zero(c); //c=0
    arb_one(d); //d=1
    arb_one(h);
    arb_mul_ui(h,h,2,prec); //h=2
    
    slong k=0,mode=0; // Tanh-Sinh = 0, Exp-Sinh = 1, Sinh-Sinh = 2
    
    if(arb_lt(bb,aa)) //若 b<a，swap bounds
    {
        arb_swap(aa, bb);
        sign=-1;
    }
    
    if( arb_is_finite(aa) && arb_is_finite(bb) )
    {
        arb_add(t_s,aa,bb,prec);
        arb_sub(t_t,bb,aa,prec);
        
        arb_mul_2exp_si(c,t_s,-1); //c = (a+b)/2
        arb_mul_2exp_si(d,t_t,-1); //d = (b-a)/2
        arb_set(v,c); //v=c
    }else if(arb_is_finite(aa))
    {
        mode=1; // Exp-Sinh
        
        //这里 ，可以续续使用前面的d=1的值
        //也可以使用 exp_sinh_opt_d 函数得到优化后的d值
        //使用优化后的d值，可以明显地加快收敛速度
        arb_set(t_s,d);
        exp_sinh_opt_d(d, func, param, order,
                       aa, eps, t_s, prec);
        
        arb_set(c,aa); //c=a
        arb_add(v,aa,d,prec); //v=a+d
    }else if(arb_is_finite(bb))
    {
        mode=1; // Exp-Sinh
        arb_neg(d,d); //d=-d
        
        //这里 ，可以继续使用前面的d=-1的值
        //也可以使用 exp_sinh_opt_d 函数得到优化后的d值
        //使用优化后的d值，可以明显地加快收敛速度
        arb_set(t_s,d);
        exp_sinh_opt_d(d, func, param, order,
                       aa, eps, t_s, prec);
        
        sign=-sign;
        arb_set(c,bb); //c=b
        arb_add(v,bb,d,prec); //v=b+d
    }else
    {
        mode=2; // Sinh-Sinh
        arb_zero(v);
    }
    
    func(s,v,param,order,prec); //s=f(v)
    
    do{
        
        //设定每次循环的初始条件
        arb_zero(p); //p=0
        arb_zero(fp); //fp=0
        arb_zero(fm); //fm=0
        
        
        arb_mul_2exp_si(h,h,-1); //h=h/2
        arb_exp(eh,h,prec); //t=eh=exp(h)
        arb_set(t,eh);
        
        if(k>0)
        {
            arb_sqr(eh,eh,prec); //eh=eh*eh
        }
        
        if(mode==0) // Tanh-Sinh 
        {
            do{
               arb_inv(t_w,t,prec); //u=exp(1/t-t)= exp(-2*sinh(j*h)) = 1/exp(sinh(j*h))^2
               arb_sub(t_s,t_w,t,prec);
               arb_exp(u,t_s,prec);
               
               arb_add_ui(t_t,u,1,prec); //r=2*u/(1+u)= 1 - tanh(sinh(j*h)) 
               arb_div(t_s,u,t_t,prec);
               arb_mul_ui(r,t_s,2,prec);
               
               
               arb_add(t_s,t_w,t,prec); //w=(t+1/t)*r/(1+u)= cosh(j*h)/cosh(sinh(j*h))^2
               arb_mul(t_s,t_s,r,prec);
               arb_div(w,t_s,t_t,prec);
               
               arb_mul(x,d,r,prec); //x=d*r
               
               arb_add(t_s,aa,x,prec);
               if( arb_gt(t_s,aa) ) //a+x > a  // if too close to a then reuse previous fp
               {
                   func(y,t_s,param,order,prec); //y=f(a+x)
                   if( arb_is_finite(y) )
                   {
                      arb_set(fp,y); // if f(x) is finite, add to local sum
                   }
               }
               
               arb_sub(t_s,bb,x,prec);
               if( arb_lt(t_s,bb) ) //b-x<b // if too close to a then reuse previous fp
               {
                   func(y,t_s,param,order,prec); //y=f(b-x)
                   if( arb_is_finite(y) )
                   {
                       arb_set(fm,y); // if f(x) is finite, add to local sum
                   }
               }
               
               arb_add(t_s,fp,fm,prec); //q = w*(fp+fm)
               arb_mul(q,t_s,w,prec);
               
               arb_add(p,p,q,prec); //p=p+q
               arb_mul(t,t,eh,prec); //t=t*eh
               
               arb_abs(t_s,q);
               arb_abs(t_t,p);
               arb_mul(t_t,t_t,eps,prec);
               
            }while( arb_gt(t_s,t_t) ); //fabs(q) > eps*fabs(p)
            
        }else
        {
            arb_mul_2exp_si(t,t,-1); //t=t/2
            
            do{
                arb_one(t_s); //t_s=1/4=0.25
                arb_div_ui(t_s,t_s,4,prec);
                
                arb_div(t_s,t_s,t,prec); //r = exp(t-.25/t) = exp(sinh(j*h))
                arb_sub(t_s,t,t_s,prec);
                arb_exp(r,t_s,prec);
                
                arb_set(w,r);//w=r
                
                arb_zero(q); //q=0
                
                if (mode==1) // Exp-Sinh
                {
                    arb_div(t_s,d,r,prec); //x = c + d/r
                    arb_add(x,t_s,c,prec);
                    
                    if(arb_eq(x,c)) //x==c // if x hit the finite endpoint then break
                    {
                        break;
                    }
                    
                    func(y,x,param,order,prec); //y=f(x)
                    
                    if(arb_is_finite(y)) // if f(x) is finite, add to local sum
                    {
                        arb_div(t_s,y,w,prec); //q += y/w
                        arb_add(q,q,t_s,prec);
                    }
                    
                }else // Sinh-Sinh
                {
                    arb_inv(t_s,r,prec); //r = (r-1/r)/2 = sinh(sinh(j*h))
                    arb_sub(t_s,r,t_s,prec);
                    arb_mul_2exp_si(r,t_s,-1);
                    
                    arb_inv(t_s,w,prec); //w = (w+1/w)/2 = cosh(sinh(j*h))
                    arb_add(t_s,t_s,w,prec);
                    arb_mul_2exp_si(w,t_s,-1);
                    
                    arb_mul(t_s,d,r,prec); //x = c - d*r
                    arb_sub(x,c,t_s,prec);
                    
                    func(y,x,param,order,prec); //y=f(x)
                    
                    if(arb_is_finite(y)) // if f(x) is finite, add to local sum
                    {
                        arb_mul(t_s,y,w,prec); //q += y*w
                        arb_add(q,q,t_s,prec);
                    }
                }
                
                arb_mul(t_s,d,r,prec); //x = c + d*r
                arb_add(x,t_s,c,prec);
                
                func(y,x,param,order,prec); //y=f(x)
                
                if(arb_is_finite(y)) // if f(x) is finite, add to local sum
                {
                    arb_mul(t_s,y,w,prec); //q += y*w
                    arb_add(q,q,t_s,prec);
                }
                
                arb_one(t_s); //t_s=1/4=0.25
                arb_div_ui(t_s,t_s,4,prec);
                
                arb_div(t_s,t_s,t,prec); //q *= t+.25/t
                arb_add(t_s,t,t_s,prec);
                arb_mul(q,q,t_s,prec);
                
                arb_add(p,p,q,prec); //p += q
                
                arb_mul(t,t,eh,prec); //t *= eh
                
                arb_abs(t_s,q);
                arb_abs(t_t,p);
                arb_mul(t_t,t_t,eps,prec);
                
            }while( arb_gt(t_s,t_t) ); //fabs(q) > eps*fabs(p)
        }
        
        arb_sub(v,s,p,prec);//v=s-p
        arb_add(s,s,p,prec); //s+=p
        ++k;
        
        arb_abs(t_s,v);
        arb_abs(t_t,s);
        arb_mul(t_t,tol,t_t,prec);
        
    } while( k<=step_min || (arb_gt(t_s,t_t) && k<=step_max) ); // k<=step_min || (fabs(v) > tol*fabs(s) && k <= step_max)
    
    /*
    arb_abs(t_s,v); //e = fabs(v)/(fabs(s)+eps) // result with estimated relative error e
    arb_abs(t_t,s);
    arb_add(t_t,t_t,eps,prec);
    arb_div(e,t_s,t_t,prec);
    */
    
    arb_mul(t_s,d,s,prec); //res=sign*d*s*h
    arb_mul(t_s,t_s,h,prec);
    arb_mul_si(res,t_s,sign,prec);
    
    
    arb_clear(tol);
    arb_clear(c);
    arb_clear(d);
    arb_clear(s);
    arb_clear(e);
    arb_clear(v);
    arb_clear(h);
    arb_clear(aa);
    arb_clear(bb);
    arb_clear(p);
    arb_clear(q);
    arb_clear(fp);
    arb_clear(fm);
    arb_clear(t);
    arb_clear(eh);
    arb_clear(u);
    arb_clear(r);
    arb_clear(w);
    arb_clear(x);
    arb_clear(y);
    arb_clear(t_s);
    arb_clear(t_t);
    arb_clear(t_w);
    
    if(k<=step_max)
        return 0;
    else
        return 1;
}



int Integration_arb(arb_t res, my_calc_func func, void *param, const slong order,
                                        const arb_t a, const arb_t b, const arb_t error,
                                        slong step_min , slong step_max,
                                        slong prec)
{
    arb_t t,gap,error_cut;
    arb_init(t);
    arb_init(gap);
    arb_init(error_cut);
    
    int res_judge=0;
    
    //当区间过大后（如功率谱过宽），在大区间上的积分很难收收敛，计算太慢
    //将给定的区间分割为多个小区间，然后并行计算
    //其中的并行计算采用 OpenMP 实现
    
    arb_sub(gap,b,a,prec);
    arb_div_ui(gap,gap,Multithreaded_divide_integration_interval_number,prec);
    
    arb_div_ui(error_cut,error,Multithreaded_divide_integration_interval_number,prec);
    
    arb_ptr sum;
    sum = _arb_vec_init(Multithreaded_divide_integration_interval_number);
    
    //Multithreaded_divide_integration_interval_number
    //Multithreaded_number
    #pragma omp parallel for firstprivate(Integral_method,step_min,step_max,prec) num_threads(Multithreaded_number) schedule(dynamic) //采用动态调度，可平衡各个线程的计算量
    for(slong i=0; i < Multithreaded_divide_integration_interval_number; i++)
    {
        arb_t a_cut,b_cut;
        arb_init(a_cut);
        arb_init(b_cut);
        
        arb_mul_ui(a_cut,gap,i,prec);
        arb_add(a_cut,a_cut,a,prec);
        
        arb_mul_ui(b_cut,gap,i+1,prec);
        arb_add(b_cut,b_cut,a,prec);
        
        
        //根据需要采用不同的积分方法
        switch(Integral_method)
        {
            case gauss_kronrod_iterate : //gauss_kronrod迭代版本
                get_gauss_kronrod_node_weight(65,prec); //获取 gauss_kronrod 的节点位置和权重
                res_judge=integration_gauss_kronrod_iterate(sum+i, func, param, order, a_cut, b_cut, error_cut, step_min, step_max, prec);
                break;
                
            case gauss_kronrod_recursive ://gauss_kronrod递归版本
                get_gauss_kronrod_node_weight(65,prec); //获取 gauss_kronrod 的节点位置和权重
                res_judge=integration_gauss_kronrod_recursive(sum+i, func, param, order, a_cut, b_cut, error_cut, step_min, step_max, prec);
                break;
                
            case double_exponential :
                res_judge=Double_Exponential_Quadrature(sum+i, func, param, order, a_cut, b_cut, error_cut, step_min, step_max, prec);
                break;
                
            default : //默认使用gauss_kronrod_iterate积分
                get_gauss_kronrod_node_weight(65,prec); //获取 gauss_kronrod 的节点位置和权重
                res_judge=integration_gauss_kronrod_iterate(sum+i, func, param, order, a_cut, b_cut, error_cut, step_min, step_max, prec);
                
        }
        {
        #pragma omp atomic //这里，需要判断积分是否达到精度要求，通过求和方式判定，全部为零则达到要求
        res_judge+=res_judge;
        }
        arb_clear(a_cut);
        arb_clear(b_cut);
    }
    
    arb_zero(t);
    for(slong i=0; i<Multithreaded_divide_integration_interval_number; i++) //求和得到积结果
    {
        arb_add(t,t,sum+i,prec);
    }
    
    arb_set(res,t);
    
    arb_clear(t);
    arb_clear(gap);
    arb_clear(error_cut);
    _arb_vec_clear(sum, Multithreaded_divide_integration_interval_number);
    
    
    if(res_judge==0)
    {
        return 0;
    }else
    {
        return 1;
    }
}



