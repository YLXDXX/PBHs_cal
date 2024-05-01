#include "find_root.h" 
#include <stdlib.h>
#include <stdbool.h>

enum FIND_ROOT_ENUM_TYPE Find_root_method; //找根方法

//Brent’s method 求根
int Find_root_brent_method(arb_t res, my_calc_func func, void *param, const slong order,
                    const arb_t a, const arb_t b, const arb_t error,
                    slong prec)
{
    //来自于 https://mathsfromnothing.au/brents-method/ 
    
    arb_t xtol,ytol,xb_1,xb_2,yb_1,yb_2,sign,c,yc,s,ys,d,x,y,t_0,t_1,t_2,t_3,t_4,t_5,t_6;
    
    // xtol,ytol 是误差项，用于控制精度
    // xb_1,xb_2 是根所在区间的端点 [xb_1,xb_2]
    // yb_1,yb_2 是根所在区间的端点的函数值 [yb_1,yb_2]
    // sign 用于符号判断 正 负 零
    // c,yc  s,ys  d 属于中间计算所用变量
    // x,y 最终结果输出 x: 根的值， y:=func(x)
    
    arb_init(xtol);
    arb_init(ytol);
    arb_init(xb_1);
    arb_init(xb_2);
    arb_init(yb_1);
    arb_init(yb_2);
    arb_init(sign);
    arb_init(c);
    arb_init(yc);
    arb_init(s);
    arb_init(ys);
    arb_init(d);
    arb_init(x);
    arb_init(y);
    arb_init(t_0);//自用变量
    arb_init(t_1);
    arb_init(t_2);
    arb_init(t_3);
    arb_init(t_4);
    arb_init(t_5);
    arb_init(t_6);
    
    //arb_set(xtol,error);
    arb_abs(xtol,error);
    arb_set(ytol,xtol);
    arb_mul_ui(ytol,ytol,100,prec); //设 ytol 比 xtol 少两个数量级
    
    arb_set(xb_1,a);
    arb_set(xb_2,b);
    
    func(yb_1,xb_1,param,order,prec);
    func(yb_2,xb_2,param,order,prec);
    
    
    arb_mul(sign,yb_1,yb_2,prec);
    if( arb_is_nonnegative(sign) )
    {
        printf("Math_Method -> find_root -> No sign change between the bounds.\nBrent's method requires the bounds to have different signs.\n");
        exit(1);
    }
    
    
    arb_abs(t_0,yb_1);
    arb_abs(t_1,yb_2);
    if( arb_lt(t_0,t_1) ) //abs(yb(1))<abs(yb(2)
    {
        arb_swap(xb_1, xb_2);
        arb_swap(yb_1, yb_2);
    }
    
    
    arb_set(c,xb_1);
    arb_set(yc,yb_1);
    
    bool usedbisec=true;
    
    //char *value_str; //用于将结果转为字符串，再将字符串转回结果 （消除多次年迭代的误差累积）
    
    while(true) //迭代求根
    {
        if( arb_ne(yb_1,yc) && arb_ne(yb_2,yc) ) // yb(1)!=yc && yb(2)!=yc
        {
            //s=xb(1)*yb(2)*yc/((yb(1)-yb(2))*(yb(1)-yc))+
            //    xb(2)*yb(1)*yc/((yb(2)-yb(1))*(yb(2)-yc))+
            //    c*yb(1)*yb(2)/((yc-yb(1))*(yc-yb(2)));
            
            arb_mul(t_0,xb_1,yb_2,prec);//xb(1)*yb(2)*yc/((yb(1)-yb(2))*(yb(1)-yc))
            arb_mul(t_0,t_0,yc,prec);
            arb_sub(t_1,yb_1,yb_2,prec);
            arb_sub(t_2,yb_1,yc,prec);
            arb_mul(t_1,t_1,t_2,prec);
            arb_div(t_3,t_0,t_1,prec);
            
            arb_mul(t_0,xb_2,yb_1,prec);//xb(2)*yb(1)*yc/((yb(2)-yb(1))*(yb(2)-yc))
            arb_mul(t_0,t_0,yc,prec);
            arb_sub(t_1,yb_2,yb_1,prec);
            arb_sub(t_2,yb_2,yc,prec);
            arb_mul(t_1,t_1,t_2,prec);
            arb_div(t_0,t_0,t_1,prec);
            arb_add(t_3,t_3,t_0,prec);
            
            arb_mul(t_0,c,yb_1,prec); //c*yb(1)*yb(2)/((yc-yb(1))*(yc-yb(2)))
            arb_mul(t_0,t_0,yb_2,prec);
            arb_sub(t_1,yc,yb_1,prec);
            arb_sub(t_2,yc,yb_2,prec);
            arb_mul(t_1,t_1,t_2,prec);
            arb_div(t_0,t_0,t_1,prec);
            
            arb_add(s,t_3,t_0,prec);
            
        }else
        {
            //s=(yb(1)*xb(2)-yb(2)*xb(1))/(yb(1)-yb(2));
            
            arb_mul(t_0,yb_1,xb_2,prec);
            arb_mul(t_1,yb_2,xb_1,prec);
            
            arb_sub(t_0,t_0,t_1,prec);
            
            
            
            arb_sub(t_1,yb_1,yb_2,prec);
            arb_div(s,t_0,t_1,prec);
            
        }
        
        
        //arb_printn(s, 40,0);printf("\n");
        
        if(!arb_is_finite(s))
        {
            //有时对于erro有过高要求，由于计算设定精度不足，会导致出现 NaN 的结果
            printf("Math_Method -> find_root -> 求根中的未知错误，可能原因：计算精度太低或函数有奇点\n可尝试使用折半法来找根\n\n");
            exit(1);
        }
        
        
        //(s-(3*xb(1)+xb(2))/4)*(s-xb(2))>=0||...
        //(usedbisec&&abs(s-xb(2))>=abs(xb(2)-c)/2)||...
        //(!usedbisec&&abs(s-xb(2))>=abs(c-d)/2)||...
        //(usedbisec&&abs(xb(2)-c)<abs(xtol))||...
        //(!usedbisec&&abs(c-d)<abs(xtol))
        
        // (s-(3*xb(1)+xb(2))/4)*(s-xb(2)) 
        arb_mul_si(t_0,xb_1,3,prec);//(s-(3*xb(1)+xb(2))/4)
        arb_add(t_0,t_0,xb_2,prec);
        arb_div_ui(t_0,t_0,4,prec);
        arb_sub(t_0,s,t_0,prec);
        
        arb_sub(t_1,s,xb_2,prec);//(s-xb(2)) 
        arb_mul(t_0,t_0,t_1,prec);
        
        // abs(s-xb(2)) 
        arb_sub(t_1,s,xb_2,prec);
        arb_abs(t_1,t_1);
        
        // abs(xb(2)-c)/2)
        arb_sub(t_2,xb_2,c,prec);
        arb_div_ui(t_2,t_2,2,prec);
        arb_abs(t_2,t_2);
        
        
        // abs(c-d)/2)
        arb_sub(t_4,c,d,prec);
        arb_div_ui(t_4,t_4,2,prec);
        arb_abs(t_4,t_4);
        
        //abs(xb(2)-c)
        arb_sub(t_5,xb_2,c,prec);
        arb_abs(t_5,t_5);
        
        //abs(c-d)
        arb_sub(t_6,c,d,prec);
        arb_abs(t_6,t_6);
        
        if( //!arb_is_finite(s) || //有时候，由于精度的问题，上一步的 s 可能会为 NaN 的结果
            arb_is_nonnegative(t_0) || //(s-(3*xb(1)+xb(2))/4)*(s-xb(2))>=0||...
            (usedbisec && arb_ge(t_1,t_2))|| //(usedbisec&&abs(s-xb(2))>=abs(xb(2)-c)/2)||...
            (!usedbisec && arb_ge(t_1,t_4)) || //(!usedbisec&&abs(s-xb(2))>=abs(c-d)/2)||...
            (usedbisec && arb_lt(t_5,xtol)) || //(usedbisec&&abs(xb(2)-c)<abs(xtol))||...
            (!usedbisec&& arb_lt(t_6,xtol)) //(!usedbisec&&abs(c-d)<abs(xtol))
        )
        {
            //s=(xb(1)+xb(2))/2;
            
            arb_add(s,xb_1,xb_2,prec);
            arb_div_ui(s,s,2,prec);
            
            usedbisec=true;
            
        }else
        {
            usedbisec=false;
        }
        
        
        //ys=f(s);
        //消除多次迭代的误差累积，先转为字符串再回来
        //value_str=arb_get_str(s,prec,ARB_STR_NO_RADIUS);
        //arb_set_str(s,value_str,prec);
        //用新方法，直接设置为中点，不累积误差
        arb_get_mid_arb(s,s); //消误差
        func(ys,s,param,order,prec);
        arb_get_mid_arb(ys,ys); //消误差
        
        arb_set(d,c); //d=c;
        arb_set(c,xb_2); //c=xb(2);
        arb_set(yc,yb_2); //yc=yb(2);
        
        
        arb_mul(sign,yb_1,ys,prec);
        if( arb_is_negative(sign) )
        {
            arb_set(xb_2,s); //xb(2)=s;
            arb_set(yb_2,ys); //yb(2)=ys;
        }else
        {
            arb_set(xb_1,s); //xb(1)=s;
            arb_set(yb_1,ys); //yb(1)=ys;
        }
        
        
        arb_abs(t_0,yb_1);
        arb_abs(t_1,yb_2);
        if( arb_lt(t_0,t_1) ) //abs(yb(1))<abs(yb(2))
        {
            arb_swap(xb_1, xb_2);
            arb_swap(yb_1, yb_2);
        }
        
        if(arb_is_zero(yb_2)) //yb(2)==0
        {
            arb_set(x,xb_2);
            arb_set(y,yb_2);
            break;
        }
        
        if(arb_is_zero(ys)) //ys==0
        {
            arb_set(x,s);
            arb_set(y,ys);
            break;
        }
        
        
        arb_sub(t_0,xb_2,xb_1,prec);
        arb_abs(t_0,t_0);
        
        if( arb_lt(t_0,xtol) ) //abs(xb(2)-xb(1))<xtol
        {
            arb_set(x,s);
            arb_set(y,ys);
            break;
        }
        
        
    }
    /*
    arb_abs(t_0,y);
    if (arb_gt(t_0,ytol)) //abs(y)>ytol
    {
        printf("warning: A pole or discontinuity may have been found.\t f(x): ");
        arb_printn(y, 20,0);printf("\n");
    }
    */
    
    //arb_set(res,x);
    arb_get_mid_arb(res,x); //消除由于迭代本身所需计算所引入的误差
    
    //释放变量
    arb_clear(xtol);
    arb_clear(ytol);
    arb_clear(xb_1);
    arb_clear(xb_2);
    arb_clear(yb_1);
    arb_clear(yb_2);
    arb_clear(sign);
    arb_clear(c);
    arb_clear(yc);
    arb_clear(s);
    arb_clear(ys);
    arb_clear(d);
    arb_clear(x);
    arb_clear(y);
    arb_clear(t_0);//自用变量
    arb_clear(t_1);
    arb_clear(t_2);
    arb_clear(t_3);
    arb_clear(t_4);
    arb_clear(t_5);
    arb_clear(t_6);
    
    return 0;
}


//二分法找根
int Find_root_bisection_method(arb_t res, my_calc_func func, void *param, const slong order,
                               const arb_t a, const arb_t b, const arb_t error,
                               slong prec)
{
    arb_t xtol,ytol,xb_1,xb_2,yb_1,yb_2,sign,s,ys,x,y,t_0,t_1;
    
    // xtol,ytol 是误差项，用于控制精度
    // xb_1,xb_2 是根所在区间的端点 [xb_1,xb_2]
    // yb_1,yb_2 是根所在区间的端点的函数值 [yb_1,yb_2]
    // sign 用于符号判断 正 负 零
    // c,yc  s,ys  d 属于中间计算所用变量
    // x,y 最终结果输出 x: 根的值， y:=func(x)
    
    arb_init(xtol);
    arb_init(ytol);
    arb_init(xb_1);
    arb_init(xb_2);
    arb_init(yb_1);
    arb_init(yb_2);
    arb_init(sign);
    arb_init(s);
    arb_init(ys);
    arb_init(x);
    arb_init(y);
    arb_init(t_0);//自用变量
    arb_init(t_1);
    
    
    //arb_set(xtol,error);
    arb_abs(xtol,error);
    arb_set(ytol,xtol);
    arb_mul_ui(ytol,ytol,100,prec); //设 ytol 比 xtol 少两个数量级
    
    arb_set(xb_1,a);
    arb_set(xb_2,b);
    
    func(yb_1,xb_1,param,order,prec);
    func(yb_2,xb_2,param,order,prec);
    
    
    arb_mul(sign,yb_1,yb_2,prec);
    if( arb_is_nonnegative(sign) )
    {
        printf("Math_Method -> find_root -> No sign change between the bounds.\nBisection's method requires the bounds to have different signs.\n");
        exit(1);
    }
    
    
    arb_abs(t_0,yb_1);
    arb_abs(t_1,yb_2);
    if( arb_lt(t_0,t_1) ) //abs(yb(1))<abs(yb(2)
    {
        arb_swap(xb_1, xb_2);
        arb_swap(yb_1, yb_2);
    }
    
    
    
    
    while( true ) //迭代求根
    {
        
        //s=(xb(1)+xb(2))/2;
        arb_add(s,xb_1,xb_2,prec);
        arb_div_ui(s,s,2,prec);
        
        //ys=f(s);
        func(ys,s,param,order,prec);
        
        
        arb_mul(sign,yb_1,ys,prec);
        if( arb_is_negative(sign) )
        {
            arb_set(xb_2,s); //xb(2)=s;
            arb_set(yb_2,ys); //yb(2)=ys;
        }else
        {
            arb_set(xb_1,s); //xb(1)=s;
            arb_set(yb_1,ys); //yb(1)=ys;
        }
        
        
        arb_abs(t_0,yb_1);
        arb_abs(t_1,yb_2);
        if( arb_lt(t_0,t_1) ) //abs(yb(1))<abs(yb(2))
        {
            arb_swap(xb_1, xb_2);
            arb_swap(yb_1, yb_2);
        }
        
        if(arb_is_zero(yb_2)) //yb(2)==0
        {
            arb_set(x,xb_2);
            arb_set(y,yb_2);
            break;
        }
        
        if(arb_is_zero(ys)) //ys==0
        {
            arb_set(x,s);
            arb_set(y,ys);
            break;
        }
        
        
        arb_sub(t_0,xb_2,xb_1,prec);
        arb_abs(t_0,t_0);
        
        
        if( arb_lt(t_0,xtol) ) //abs(xb(2)-xb(1))<xtol
        {
            arb_set(x,s);
            arb_set(y,ys);
            break;
        }
        
        
    }
    
    /*
    arb_abs(t_0,y);
    if (arb_gt(t_0,ytol)) //abs(y)>ytol
    {
        printf("warning: A pole or discontinuity may have been found.\t f(x): ");
        arb_printn(y, 20,0);printf("\n");
    }
    */
    arb_set(res,x);
    
    
    //释放变量
    arb_clear(xtol);
    arb_clear(ytol);
    arb_clear(xb_1);
    arb_clear(xb_2);
    arb_clear(yb_1);
    arb_clear(yb_2);
    arb_clear(sign);
    arb_clear(s);
    arb_clear(ys);
    arb_clear(x);
    arb_clear(y);
    arb_clear(t_0);//自用变量
    arb_clear(t_1);
    
    return 0;
}


//在指定的区间内找根算法
//指定函数，指定区间，指定查找子区间个数，指定精度
int Find_interval_root(arb_t res, my_calc_func func, void *param, const slong order,
                       const arb_t a, const arb_t b, const arb_t error,
                       const slong num, const enum Root_type r_type, slong prec)
{
    //由于这里，可能有多个根，需要求最小的根
    //将区间 [a,b] 分为 num 份小区间，依次查看是否有根
    //当知道某个小区间有根后，再用求根算法求解
    
    //对于找C(r)的最大值，需要特殊处理，对于有些特殊的profile，第一个极值点可能是最小值点
    arb_ptr find_point; //数组，用于存储各个区间点
    find_point=_arb_vec_init(num+1); //分配相应的空间
    
    
    //局部变量
    arb_t t,aa,bb,mid,gap,front_refer,after_refer;
    
    arb_init(t);
    arb_init(aa);
    arb_init(bb);
    arb_init(mid);
    arb_init(gap);
    arb_init(front_refer);
    arb_init(after_refer);
    
    
    arb_set(aa,a);
    arb_set(bb,b);
    
    arb_sub(gap,bb,aa,prec);
    arb_div_si(gap,gap,num,prec); //gap=(b-a)/num
    
    func(front_refer,aa,param,order,prec);
    
    arb_set(find_point,front_refer);//备用
    
    
    int judge=0;
    
    for(slong i=1; i<=num ; i++ )
    {
        //当前的区间位置 a+i*gap
        arb_mul_si(mid,gap,i,prec);
        arb_add(mid,mid,a,prec);
        
        func(after_refer,mid,param,order,prec);
        
        arb_set(find_point+i,after_refer); //备用
        
        //不用两个数相乘为负得到有根，有时候会有不在定义域值被设为零的情况
        if( (arb_is_negative(find_point+i-1) && arb_is_positive(find_point+i)) ||
            (arb_is_positive(find_point+i-1) && arb_is_negative(find_point+i))
          )
        {
            //找C(r)的最大值，利用一阶导数，前一个值为正，后一个值为负
            if(r_type==Root_C_Max)
            {
                //注意，ζ’ + rζ’’ 前面有个负号，正负号刚好相反
                if(arb_is_negative(find_point+i-1) && arb_is_positive(find_point+i) )
                {
                   //找到C(r)最大值的位置 r_m，结束查找
                    arb_set(bb,mid); //先设置bb
                    
                    arb_mul_si(mid,gap,i-1,prec);
                    arb_add(aa,mid,a,prec);
                    
                    judge=1;
                    
                    break;
                }else
                {
                    //找到的是C(r)最小值的位置，继续向后查找
                    //arb_set(front_refer,find_point+i);
                    
                    arb_set(t,mid);
                    arb_mul_si(mid,gap,i-1,prec);
                    arb_add(mid,mid,a,prec);
                    
                    if(Stdout_verbose==true)
                    {
                        printf("\n\n\t 发现最小值所在区间：[");
                        arb_printn(mid, 15,0);printf(",");
                        arb_printn(t, 15,0);printf("]\n");
                        printf("\t 区间端处的值为：[");
                        arb_printn(find_point+i-1, 15,0);printf(",");
                        arb_printn(find_point+i, 15,0);printf("]\n\n");
                    }
                    
                    continue;
                }
            }
            
            
            arb_set(bb,mid); //先设置bb
            
            arb_mul_si(mid,gap,i-1,prec);
            arb_add(aa,mid,a,prec);
            
            judge=1;
            
            break;
        }
        
    }
    
    //arb_printn(aa, 15,0);printf("\n");
    //arb_printn(bb, 15,0);printf("\n\n");
    if( judge==1 )
    {
        //利用求根算法求根，有两种方法可选，建议使用 brent's method
        //根据需要采用不同的积分方法
        switch(Find_root_method)
        {
            case Brent_method :
                Find_root_brent_method(t, func, param, order,aa, bb, error, prec);
                break;
            case bisection_methode :
                Find_root_bisection_method(t, func, param, order,aa, bb, error, prec);
                break;
            default :
                printf("Math_Method -> find_root -> Find_root_method 输入有误\n");
                exit(1);
        }
        
        arb_set(res,t);
        
    }else
    {
        printf("Math_Method -> find_root -> 在所给区间和所给step内，未找到根\n[");
        arb_printn(a, 7,0);printf(","); //打印变量
        arb_printn(b, 7,0);printf("]\tstep: "); //打印变量
        printf("%li\n",num); //打印变量
        exit(1); 
    }
    
    
    _arb_vec_clear(find_point, num+1);
    arb_clear(t);
    arb_clear(aa);
    arb_clear(bb);
    arb_clear(mid);
    arb_clear(gap);
    arb_clear(front_refer);
    arb_clear(after_refer);
    
    return 0;
}




//在一个给定的范围内，尽可能的找出该范围内的所有根
int Find_interval_multi_root(arb_ptr* res, my_calc_func func, void *param, const slong order,
                       const arb_t a, const arb_t b, const arb_t error,
                       const slong num, slong prec)
{
    //区间 [a,b] 可能有多个根，尽可能需求所有根
    //将 [a,b] 分为 num 份小区间，依次检查每个小区间是否有根
    //当知道某个小区间有根后，再用求根算法求解
    //返回结果是根的数目
    
    //注意，这里的结果以数组形式返回，数组 res 在传入的时候不要初始化，后面会统一初始化
    //调用示例
    //arb_ptr muil_r; //muil_r是一个指针，
    //arb_ptr* m_r; //改变muil_r指针指向的地址，需要一个指向该指针的指针
    //m_r=&muil_r;
    //将 m_r 作为结果传入函数即可
    
    arb_ptr find_point; //数组，用于存储各个区间端点
    find_point=_arb_vec_init(num+1); //分配相应的空间
    
    
    arb_ptr multi_root; //数组，用于存储多个根
    multi_root=_arb_vec_init(1); //由于事先不知道有多少个根，先分配1个，后面不够再增加
    
    slong root_num; //根计数
    root_num=0;
    
    //局部变量
    arb_t t,aa,bb,mid,gap;
    
    arb_init(t);
    arb_init(aa);
    arb_init(bb);
    arb_init(mid);
    arb_init(gap);
    
    
    arb_set(aa,a);
    arb_set(bb,b);
    
    arb_sub(gap,bb,aa,prec);
    arb_div_si(gap,gap,num,prec); //gap=(b-a)/num
    
    func(t,aa,param,order,prec);
    arb_set(find_point,t);//第一个点，编号0
    
    
    for(slong i=1; i<=num ; i++ )
    {
        //当前的区间位置 a+i*gap
        arb_mul_si(mid,gap,i,prec);
        arb_add(mid,mid,a,prec);
        
        func(t,mid,param,order,prec);
        arb_set(find_point+i,t); //第i+1个点
        
        
        //判断 [i-1,i]这个小区间内是否有根
        if( (arb_is_negative(find_point+i-1) && arb_is_positive(find_point+i)) ||
            (arb_is_positive(find_point+i-1) && arb_is_negative(find_point+i))
        )
        {
            //有根，找根，然后存入 multi_root 
            //区间[aa,bb]
            arb_set(bb,mid); //先设置bb
            
            arb_mul_si(mid,gap,i-1,prec);
            arb_add(aa,mid,a,prec);
            
            //利用求根算法求根，有两种方法可选，建议使用 brent's method
            switch(Find_root_method)
            {
                case Brent_method :
                    Find_root_brent_method(t, func, param, order,aa, bb, error, prec);
                    break;
                case bisection_methode :
                    Find_root_bisection_method(t, func, param, order,aa, bb, error, prec);
                    break;
                default :
                    printf("Math_Method -> find_root -> Find_root_method 输入有误\n");
                    exit(1);
            }
            
            arb_set(multi_root+root_num,t);//第root_num+1个根，root_num初始值为0
            
            root_num=root_num+1;
            
            //为下一个根分配空间
            multi_root=flint_realloc(multi_root, (root_num+1) * sizeof(arb_struct));
            
            arb_init(multi_root+root_num); //新分配的数组空间，需要手动初始化
            
        }
        
    }
    
    //arb_printn(aa, 15,0);printf("\n");
    //arb_printn(bb, 15,0);printf("\n\n");
    if( root_num==0 )
    {
        printf("Math_Method -> find_root -> 在所给区间和所给step内，未找到根\n[");
        arb_printn(a, 7,0);printf(","); //打印变量
        arb_printn(b, 7,0);printf("]\tstep: "); //打印变量
        printf("%li\n",num); //打印变量
        exit(1); 
        
    }else
    {
        //找到的根以数组的形式返回，传出的是数组的地址
        arb_ptr t_res;
        t_res=NULL;
        t_res=flint_realloc(t_res, root_num*sizeof(arb_struct)); //手动初始化，根的个数为 root_num
        //初始化 res 后赋值
        for (slong i=0; i < root_num; i++)
        {
            arb_init(t_res+i); //分配空间后初始化
            arb_set(t_res+i, multi_root+i); //赋值
            
        }
        
        *res=t_res; //将t_res的地址传出，在本函数内不清理t_res
        
    }
    
    
    _arb_vec_clear(find_point, num+1); //清理数组
    _arb_vec_clear(multi_root, root_num+1); //清理数组
    arb_clear(t);
    arb_clear(aa);
    arb_clear(bb);
    arb_clear(mid);
    arb_clear(gap);
    
    return root_num; //返回根的个数
}

