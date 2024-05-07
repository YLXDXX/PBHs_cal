#include "quadrature_binary_rectangle_double_exponential.h"
#include "quadrature.h"
#include <stdlib.h>

// return optimized Exp-Sinh integral split point d
static void exp_sinh_opt_d(arb_t res, const arb_t f_x,
                           my_calc_func_binary func, void *param, const slong order,
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
    func(t_w,f_x,t_t,param,order,prec); 
    arb_mul_ui(t_w,t_w,4,prec);
    
    arb_mul_2exp_si(t_s,d,-1); //f(a + d/2)
    arb_add(t_s,t_s,a,prec);
    func(t_t,f_x,t_s,param,order,prec);
    
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
            func(fl,f_x,t_s,param,order,prec);
            
            arb_mul(t_s,d,r,prec); //fr = f(a + d*r)*r*r
            arb_add(t_s,t_s,a,prec);
            func(fr,f_x,t_s,param,order,prec);
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
                func(fl,f_x,t_s,param,order,prec);
                
                arb_mul(t_s,d,r,prec); //fr = f(a + d*r)*r*r
                arb_add(t_s,t_s,a,prec);
                func(fr,f_x,t_s,param,order,prec);
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

//二元函数积分，先对y进行积分 ∫f(x,y)dy
static int binary_integration_double_exponential_y(arb_t res, const arb_t f_x,
                                            my_calc_func_binary func, void *param, const slong order,
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
        exp_sinh_opt_d(d, f_x, func, param, order,
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
        exp_sinh_opt_d(d, f_x, func, param, order,
                       aa, eps, t_s, prec);
        
        sign=-sign;
        arb_set(c,bb); //c=b
        arb_add(v,bb,d,prec); //v=b+d
    }else
    {
        mode=2; // Sinh-Sinh
        arb_zero(v);
    }
    
    func(s,f_x,v,param,order,prec); //s=f(v)
    
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
                    func(y,f_x,t_s,param,order,prec); //y=f(a+x)
                    if( arb_is_finite(y) )
                    {
                        arb_set(fp,y); // if f(x) is finite, add to local sum
                    }
                }
                
                arb_sub(t_s,bb,x,prec);
                if( arb_lt(t_s,bb) ) //b-x<b // if too close to a then reuse previous fp
                {
                    func(y,f_x,t_s,param,order,prec); //y=f(b-x)
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
                    
                    func(y,f_x,x,param,order,prec); //y=f(x)
                    
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
                    
                    func(y,f_x,x,param,order,prec); //y=f(x)
                    
                    if(arb_is_finite(y)) // if f(x) is finite, add to local sum
                    {
                        arb_mul(t_s,y,w,prec); //q += y*w
                        arb_add(q,q,t_s,prec);
                    }
                }
                
                arb_mul(t_s,d,r,prec); //x = c + d*r
                arb_add(x,t_s,c,prec);
                
                func(y,f_x,x,param,order,prec); //y=f(x)
                
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
     *    arb_abs(t_s,v); //e = fabs(v)/(fabs(s)+eps) // result with estimated relative error e
     *    arb_abs(t_t,s);
     *    arb_add(t_t,t_t,eps,prec);
     *    arb_div(e,t_s,t_t,prec);
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




// return optimized Exp-Sinh integral split point d
static void exp_sinh_opt_d_integral_x(arb_t res,
                                      my_calc_func_binary func, void *param, const slong order,
                                      const arb_t y_a, const arb_t y_b, const arb_t y_error,
                                      slong y_step_min , slong y_step_max,
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
    //func(t_w,t_t,param,order,prec); 
    binary_integration_double_exponential_y(t_w,t_t,func,param,order,
                                               y_a,y_b,y_error,y_step_min,y_step_max,prec);
    
    arb_mul_ui(t_w,t_w,4,prec);
    
    arb_mul_2exp_si(t_s,d,-1); //f(a + d/2)
    arb_add(t_s,t_s,a,prec);
    //func(t_t,t_s,param,order,prec);
    binary_integration_double_exponential_y(t_t,t_s,func,param,order,
                                               y_a,y_b,y_error,y_step_min,y_step_max,prec);
    
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
            //func(fl,t_s,param,order,prec);
            binary_integration_double_exponential_y(fl,t_s,func,param,order,
                                                       y_a,y_b,y_error,y_step_min,y_step_max,prec);
            
            arb_mul(t_s,d,r,prec); //fr = f(a + d*r)*r*r
            arb_add(t_s,t_s,a,prec);
            //func(fr,t_s,param,order,prec);
            binary_integration_double_exponential_y(fr,t_s,func,param,order,
                                                       y_a,y_b,y_error,y_step_min,y_step_max,prec);
            
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
                //func(fl,t_s,param,order,prec);
                binary_integration_double_exponential_y(fl,t_s,func,param,order,
                                                           y_a,y_b,y_error,y_step_min,y_step_max,prec);
                
                
                arb_mul(t_s,d,r,prec); //fr = f(a + d*r)*r*r
                arb_add(t_s,t_s,a,prec);
                //func(fr,t_s,param,order,prec);
                binary_integration_double_exponential_y(fr,t_s,func,param,order,
                                                           y_a,y_b,y_error,y_step_min,y_step_max,prec);
                
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



//对y积分后，再对x积分 ∫[∫f(x,y)dy]dx
static int binary_integration_double_exponential_y_x(arb_t res,
                                              my_calc_func_binary func, void *param, const slong order,
                                              const arb_t y_a, const arb_t y_b, const arb_t y_error,
                                              slong y_step_min , slong y_step_max,
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
        exp_sinh_opt_d_integral_x(d, func, param, order,
                                  y_a, y_b, y_error,
                                  y_step_min, y_step_max,
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
        exp_sinh_opt_d_integral_x(d, func, param, order,
                                  y_a, y_b, y_error,
                                  y_step_min, y_step_max,
                                  aa, eps, t_s, prec);
        
        sign=-sign;
        arb_set(c,bb); //c=b
        arb_add(v,bb,d,prec); //v=b+d
    }else
    {
        mode=2; // Sinh-Sinh
        arb_zero(v);
    }
   
    //func(s,v,param,order,prec); //s=f(v)
    binary_integration_double_exponential_y(s,v,func,param,order,
                                               y_a,y_b,y_error,y_step_min,y_step_max,prec);
    
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
                    //func(y,t_s,param,order,prec); //y=f(a+x)
                    binary_integration_double_exponential_y(y,t_s,func,param,order,
                                                               y_a,y_b,y_error,y_step_min,y_step_max,prec);
                    
                    if( arb_is_finite(y) )
                    {
                        arb_set(fp,y); // if f(x) is finite, add to local sum
                    }
                }
                
                arb_sub(t_s,bb,x,prec);
                if( arb_lt(t_s,bb) ) //b-x<b // if too close to a then reuse previous fp
                {
                    //func(y,t_s,param,order,prec); //y=f(b-x)
                    binary_integration_double_exponential_y(y,t_s,func,param,order,
                                                               y_a,y_b,y_error,y_step_min,y_step_max,prec);
                    
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
                    
                    //func(y,x,param,order,prec); //y=f(x)
                    binary_integration_double_exponential_y(y,x,func,param,order,
                                                               y_a,y_b,y_error,y_step_min,y_step_max,prec);
                    
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
                    
                    //func(y,x,param,order,prec); //y=f(x)
                    binary_integration_double_exponential_y(y,x,func,param,order,
                                                               y_a,y_b,y_error,y_step_min,y_step_max,prec);
                    
                    if(arb_is_finite(y)) // if f(x) is finite, add to local sum
                    {
                        arb_mul(t_s,y,w,prec); //q += y*w
                        arb_add(q,q,t_s,prec);
                    }
                }
                
                arb_mul(t_s,d,r,prec); //x = c + d*r
                arb_add(x,t_s,c,prec);
                
                //func(y,x,param,order,prec); //y=f(x)
                binary_integration_double_exponential_y(y,x,func,param,order,
                                                           y_a,y_b,y_error,y_step_min,y_step_max,prec);
                
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
     *    arb_abs(t_s,v); //e = fabs(v)/(fabs(s)+eps) // result with estimated relative error e
     *    arb_abs(t_t,s);
     *    arb_add(t_t,t_t,eps,prec);
     *    arb_div(e,t_s,t_t,prec);
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

//二元函数积分：先对y进再对x ∫[∫f(x,y)dy]dx
int quadrature_binary_rectangle_double_exponential(arb_t res, my_calc_func_binary func, void *param, const slong order,
                                                   const arb_t x_a, const arb_t x_b, const arb_t x_error, 
                                                   slong x_step_min , slong x_step_max,
                                                   const arb_t y_a, const arb_t y_b, const arb_t y_error,
                                                   slong y_step_min , slong y_step_max,
                                                   slong prec)
{
    arb_t s;
    arb_init(s);
    
    int i;
    
    i=binary_integration_double_exponential_y_x(s,func,param,order,
                                                y_a,y_b,y_error,y_step_min,y_step_max,
                                                x_a,x_b,x_error,x_step_min,x_step_max,
                                                prec);
    
    arb_set(res,s);
    arb_clear(s);
    
    return i;
}



