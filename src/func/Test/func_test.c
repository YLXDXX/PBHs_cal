#include "func_test.h" 
#include <stdlib.h>

void Cal_test(slong prec)
{
    //输出拟合函数的点
    
    FILE *fdp ;
    
    fdp = fopen(Out_fitted_file,"w");
    if( fdp == NULL ) { //对文件打开操作进行判断
        printf("\n\nOpen Error: %s\t\n",Out_fitted_file);perror("file");printf("\n");
        exit(-1);
    }
    
    
    double a=-1;
    double b=1;
    int n=90;
    double gap=(b-a)/n;
    arb_t ss,tt;
    arb_init(ss);
    arb_init(tt);
    
    for(int i=0; i<n ; i+=1)
    {
        double gg=a+i*gap;
        
        arb_set_d(tt,gg);
        
        Fit_get_value(ss, tt, Fit_test, 0, prec);
        
        arb_fprintn(fdp,tt,8,ARB_STR_NO_RADIUS );
        fprintf(fdp,"\t");
        arb_fprintn(fdp,ss,8,ARB_STR_NO_RADIUS );
        fprintf(fdp,"\n");
    }
    
    fclose(fdp);
}



int Func_test(arb_t res, const arb_t x, void* params, const slong order, slong prec) //测试函数，用于各种测试
{
    //arb_log(res,x,prec);
    
    arb_exp(res,x,prec);
    arb_add(res,res,x,prec);
    arb_sin(res,res,prec);
    
    //arb_sqr(res,x,prec);
    //arb_mul_si(res,res,25,prec);
    //arb_add_si(res,res,1,prec);
    //arb_inv(res,res,prec);
   
    return 0;
}


int Find_root_test(arb_t res, const arb_t x, void* params, const slong order, slong prec)
{
    
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    
    //x^3-13*x^2+20*x+100
    /*
    arb_pow_ui(s,x,3,prec);
    arb_sqr(t,x,prec);
    arb_mul_ui(t,t,13,prec);
    arb_sub(s,s,t,prec);
    arb_mul_ui(t,x,20,prec);
    arb_add(s,s,t,prec);
    arb_add_ui(res,s,100,prec);
    */
    
    //x^2-2
    //arb_sqr(t,x,prec);
    //arb_sub_ui(res,t,2,prec);
    
    //ln(x)
    arb_log(res,x,prec);
    
    arb_clear(s);
    arb_clear(t);
    
    return 0;
}



//测试函数用 Func_test_3
void PR(arb_t res, const arb_t k, const arb_t A, const arb_t sigma, const arb_t kstar, slong prec)
{
    
    arb_t s,t,w;
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //右边分子
    arb_div(s,k,kstar,prec);
    arb_log(s,s,prec);
    arb_sqr(s,s,prec);
    
    arb_sqr(t,sigma,prec);
    arb_mul_ui(t,t,2,prec);
    arb_div(s,s,t,prec);
    arb_neg(s,s);
    arb_exp(s,s,prec);
    
    //右边分母
    arb_const_pi(t,prec);
    arb_mul_ui(t,t,2,prec);
    arb_sqrt(t,t,prec);
    arb_mul(t,t,sigma,prec);
    
    arb_div(s,s,t,prec);
    
    //左边系数
    arb_mul(res,s,A,prec);
    
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
}

int Func_test_2(arb_t res, const arb_t x, const arb_t y, void* params, const slong order, slong prec) //测试函数，用于各种测试
{
    //arb_one(res);
    
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    arb_add(s,x,y,prec);
    arb_neg(s,s);
    arb_exp(res,s,prec);
    
    arb_clear(s);
    arb_clear(t);
    return 0;
}

int Func_test_3(arb_t res, const arb_t p1, const arb_t p2, void* params, const slong order, slong prec) //测试函数，用于各种测试
{
    //arb_one(res);
    
    arb_t s,t,w,kR,sigma,sigmatau,kstar,A;
    
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    arb_init(kR);
    arb_init(sigma);
    arb_init(sigmatau);
    arb_init(kstar);
    arb_init(A);
    
    //基本变量赋值
    arb_one(kR);
    arb_one(sigma);
    arb_one(sigmatau);
    arb_one(kstar);
    arb_one(A);
    
    //arb_mul_si(kR,kR,100000000000000000000,prec); //10^20
    arb_set_str(kR,"1E20",prec);
    
    arb_div_ui(sigma,sigma,100,prec); //0.02
    arb_mul_ui(sigma,sigma,2,prec); 
    
    arb_div_ui(sigmatau,sigmatau,10,prec); //0.1
    
    //arb_mul_si(kstar,kstar,100000000000000000000,prec); //10^20
    arb_set(kstar,kR);
    
    arb_div_ui(A,A,10,prec); //0.1
    
    
    
    //公式从左往右，每个乘号，一个分割
    
    arb_mul(s,p1,p2,prec);
    arb_mul_ui(s,s,4,prec);
    arb_inv(s,s,prec);
    
    arb_div(t,p1,kR,prec);
    arb_sqr(w,t,prec); //后面重复用
    arb_sqr(t,w,prec);
    arb_mul(s,s,t,prec);
    
    
    arb_neg(w,w);//复用完成
    arb_exp(w,w,prec);
    arb_mul(s,s,w,prec);
    
    
    arb_div(t,p2,kR,prec);
    arb_sqr(w,t,prec); //后面重复用
    arb_sqr(t,w,prec);
    arb_mul(s,s,t,prec);
    
    
    arb_neg(w,w);//复用完成
    arb_exp(w,w,prec);
    arb_mul(s,s,w,prec);
    
    PR(t,p1,A,sigma,kstar,prec);
    arb_mul(s,s,t,prec);
    
    
    //分子
    arb_add(t,p1,p2,prec);
    arb_mul_ui(w,kstar,2,prec);
    arb_div(t,t,w,prec);
    arb_log(t,t,prec);
    arb_sqr(t,t,prec);
    
    arb_sqr(w,sigmatau,prec);
    arb_mul_ui(w,w,2,prec);
    arb_div(t,t,w,prec);
    arb_neg(t,t);
    arb_exp(t,t,prec);
    
    //分母
    arb_const_pi(w,prec);
    arb_mul_ui(w,w,2,prec);
    arb_sqrt(w,w,prec);
    arb_mul(w,w,sigmatau,prec);
    
    arb_div(t,t,w,prec);
    arb_mul(s,s,t,prec);
    
    PR(t,p2,A,sigma,kstar,prec);
    arb_mul(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    
    arb_clear(kR);
    arb_clear(sigma);
    arb_clear(sigmatau);
    arb_clear(kstar);
    arb_clear(A);
    
    return 0;
}


//一元函数多个根查找测试
int Func_test_4(arb_t res, const arb_t x, void* params, const slong order, slong prec) //测试函数，用于各种测试
{
    arb_t s,t,w;
    
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    /*
    //y=x^3-8*x^2+20 有3个根[-2,10]
    arb_pow_ui(t,x,3,prec);
    arb_pow_ui(s,x,2,prec);
    arb_mul_ui(s,s,8,prec);
    arb_sub(t,t,s,prec);
    arb_add_si(res,t,20,prec);
    */
    //y=sin(4x) 有15个根[-2,10]
    arb_mul_ui(t,x,4,prec);
    arb_sin(res,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    
    return 0; 
}


int Func_test_5(arb_t res, const arb_t x, void* params, const slong order, slong prec) //测试函数，用于各种测试
{
    arb_t s,t,w;
    
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //f(x)=x^(-0.8) 在[0,1]上积分精确值为5
    arb_set_str(s,"-0.8",prec);
    arb_pow(res,x,s,prec);
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    
    return 0; 
}

int Func_test_6(arb_t res, const arb_t x, void* params, const slong order, slong prec) //测试函数，用于各种测试
{
    arb_t s,t,w;
    
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //f(x)=e^(-0.2*x) 在[0,+∞]上积分精确值为5
    arb_set_str(s,"-0.2",prec);
    arb_mul(s,s,x,prec);
    arb_exp(res,s,prec);
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    
    return 0; 
}


int Func_test_7(arb_t res, const arb_t x, void* params, const slong order, slong prec) //测试函数，用于各种测试
{
    arb_t s,t,w;
    
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //f(x)=e^(-x^2) 
    arb_sqr(s,x,prec);
    arb_neg(s,s);
    arb_exp(res,s,prec);
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    
    return 0; 
}

int Func_test_8(arb_t res, const arb_t x, void* params, const slong order, slong prec) //测试函数，用于各种测试
{
    arb_t s,t,w;
    
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //f(x)=[e^(-0.2/x)]/y^2 在[0,+∞]上积分精确值为5
    arb_set_str(s,"-0.2",prec);
    arb_div(s,s,x,prec);
    arb_exp(s,s,prec);
    
    arb_sqr(t,x,prec);
    
    arb_div(res,s,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    
    return 0;
}

int Func_test_9(arb_t res, const arb_t x, void* params, const slong order, slong prec) //测试函数，用于各种测试
{
    arb_t s,t,w;
    
    arb_init(s);
    arb_init(t);
    arb_init(w);
    
    //振荡函数积分
    arb_set_str(t,"0.00001",prec);
    arb_mul(s,x,t,prec);
    arb_sin(s,s,prec);
    arb_sqr(s,s,prec);
    
    arb_mul(res,s,x,prec);
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(w);
    
    return 0;
}

int Func_test_10(arb_t res, const arb_t x, void* params, const slong order, slong prec) //测试函数，用于各种测试
{
    arb_t s,t;
    arb_init(s);
    arb_init(t);
    
    //发散函数积分 1/(1-x)
    arb_neg(t,x);
    arb_add_ui(t,t,1,prec);
    arb_inv(res,t,prec);
    
    arb_clear(s);
    arb_clear(t);
    return 0;
}
