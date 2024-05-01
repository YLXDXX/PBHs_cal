#include "non_gaussianity.h"
#include <stdlib.h>

// exponential_tail 类型函数及其导数
int Non_Gaussianity_exponential_tail_n(arb_t res, const arb_t x, const slong order, slong prec)
{
    arb_t s,t,beta_inverse;
    
    arb_init(s);
    arb_init(t);
    arb_init(beta_inverse);
    
    // exponential_tail 非高斯性的函数为 ζ=-1/β*ln(1-β*ζ_G)
    //即 y=-1/β*ln(1-β*x)
    
    //判断β的正负，确定 x 取值范围
    arb_inv(beta_inverse,Exponential_tail_beta,prec);
    if( arb_is_positive(Exponential_tail_beta) ) //β>0
    {
        if( arb_ge(x,beta_inverse) ) //定义域 x<1/β
        {
            arb_zero(res);
            
            arb_clear(s);
            arb_clear(t);
            arb_clear(beta_inverse);
            return 1;
        }
    }else //β<0
    {
        if( arb_le(x,beta_inverse) ) //定义域 x>1/β
        {
            arb_zero(res);
            
            arb_clear(s);
            arb_clear(t);
            arb_clear(beta_inverse);
            return 1;
        }
    }
    
    switch(order) 
    {
        case 0 : //原函数
            //(-1/β)*ln(1-β*x)
            
            arb_mul(s,x,Exponential_tail_beta,prec); // ln(1-β*x)
            arb_neg(s,s);
            arb_add_si(s,s,1,prec);
            arb_abs(s,s); //取绝对值，防止取对数遇到问题
            arb_log(s,s,prec); 
            
            arb_mul(s,s,beta_inverse,prec);
            arb_neg(res,s); //最后取负号
            
            break;
            
        case 1 : // 一阶导
            // 1/(1-β*x)
            
            arb_mul(s,x,Exponential_tail_beta,prec); //分母
            arb_neg(s,s);
            arb_add_si(s,s,1,prec);
            
            arb_inv(res,s,prec);
            
            break;
            
        case 2: // 二阶导
            //β/(1-β*x)^2
            arb_mul(s,x,Exponential_tail_beta,prec); //分母
            arb_neg(s,s);
            arb_add_si(s,s,1,prec);
            arb_sqr(s,s,prec);
            
            arb_div(res,Exponential_tail_beta,s,prec);
            
            break;
            
        case 3 : // 三阶导
            //2*(β^2)/(1-β*x)^3
            arb_mul(s,x,Exponential_tail_beta,prec); //分母
            arb_neg(s,s);
            arb_add_si(s,s,1,prec);
            arb_pow_ui(s,s,3,prec);
            
            arb_sqr(t,Exponential_tail_beta,prec);//分子
            arb_mul_si(t,t,2,prec);
            
            arb_div(res,t,s,prec);
            
            break;
            
        case 4 : // 四阶导
            // 3*2*(β^3)/(1-β*x)^4
            arb_mul(s,x,Exponential_tail_beta,prec); //分母
            arb_neg(s,s);
            arb_add_si(s,s,1,prec);
            arb_pow_ui(s,s,4,prec);
            
            arb_pow_ui(t,Exponential_tail_beta,3,prec);//分子
            arb_mul_si(t,t,6,prec);
            
            arb_div(res,t,s,prec);
            
            break;
            
        default:
            printf("General -> non_gaussianity -> Non_Gaussianity_exponential_tail_n 输入有误\n");
            exit(1);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(beta_inverse);
    
    return 0;
}


// up_step 类型函数及其导数
int Non_Gaussianity_up_step_n(arb_t res, const arb_t x, const slong order, slong prec)
{
    arb_t s,t,h,hx,sqrt_1hx,j1h;
    
    arb_init(s);
    arb_init(t);
    arb_init(h);
    arb_init(hx);
    arb_init(sqrt_1hx);
    arb_init(j1h);
    
    // up_step 非高斯性的函数为 R=-2/|h|*[sqrt(1-|h|R_G)-1]
    //即 y=-2/h*[sqrt(1-h*x)-1]
    // |Up_step_h| --> h
    arb_abs(h,Up_step_h); //h=|Up_step_h|
    
    
    //x的取值范围是 x<1/h
    arb_inv(j1h,h,prec);
    if(arb_gt(x,j1h))
    {
        arb_zero(res);
        
        arb_clear(s);
        arb_clear(t);
        arb_clear(h);
        arb_clear(hx);
        arb_clear(sqrt_1hx);
        arb_clear(j1h);
        return 1;
    }
    
    
    arb_mul(hx,h,x,prec); // h*x
    
    arb_neg(sqrt_1hx,hx);
    arb_add_si(sqrt_1hx,sqrt_1hx,1,prec);
    arb_abs(sqrt_1hx,sqrt_1hx);//取绝对值，防止开根号遇到问题
    arb_sqrt(sqrt_1hx,sqrt_1hx,prec); // sqrt(1-h*x)
    
    
    switch(order) 
    {
        case 0 : //原函数
            // (-2/h)*(sqrt(1-h*x)-1)
            
            arb_sub_si(s,sqrt_1hx,1,prec);
            arb_mul_si(s,s,2,prec);
            arb_div(s,s,h,prec);
            arb_neg(res,s); //最后取负号
            
            break;
            
        case 1 : // 一阶导
            // 1/sqrt(1-h*x)
            
            arb_inv(res,sqrt_1hx,prec);
            
            break;
            
        case 2: // 二阶导
            // -h/(sqrt(1-h*x)*(2*h*x-2))
            
            arb_mul_si(s,hx,2,prec);//分母
            arb_sub_si(s,s,2,prec);
            arb_mul(s,s,sqrt_1hx,prec);
            
            arb_div(s,h,s,prec);
            arb_neg(res,s); //最后取负号
            
            break;
            
        case 3 : // 三阶导
            // -(3*h^2*sqrt(1-h*x))/(4*h^3*x^3-12*h^2*x^2+12*h*x-4)
            
            
            arb_pow_ui(t,hx,3,prec);//分子
            arb_mul_si(t,t,4,prec);
            
            arb_sqr(s,hx,prec);
            arb_mul_si(s,s,12,prec);
            arb_sub(t,t,s,prec);
            
            arb_mul_si(s,hx,12,prec);
            arb_add(t,t,s,prec);
            
            arb_sub_si(t,t,4,prec);
            
            arb_sqr(s,h,prec); //分母
            arb_mul_si(s,s,3,prec);
            arb_mul(s,s,sqrt_1hx,prec);
            
            arb_div(s,s,t,prec);
            arb_neg(res,s); //最后取负号
            
            break;
            
        case 4 : // 四阶导
            // -(15*h^3)/(sqrt(1-h*x)*(8*h^3*x^3-24*h^2*x^2+24*h*x-8))
            
            arb_pow_ui(t,hx,3,prec);//分子
            arb_mul_si(t,t,8,prec);
            
            arb_sqr(s,hx,prec);
            arb_mul_si(s,s,24,prec);
            arb_sub(t,t,s,prec);
            
            arb_mul_si(s,hx,24,prec);
            arb_add(t,t,s,prec);
            
            arb_sub_si(t,t,8,prec);
            
            arb_mul(t,t,sqrt_1hx,prec);
            
            arb_pow_ui(s,h,3,prec);//分母
            arb_mul_si(s,s,15,prec);
            
            arb_div(s,s,t,prec);
            arb_neg(res,s); //最后取负号
            
            break;
            
        default:
            printf("General -> non_gaussianity -> Non_Gaussianity_up_step_n 输入有误\n");
            exit(1);
    }
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(h);
    arb_clear(hx);
    arb_clear(sqrt_1hx);
    arb_clear(j1h);
    
    return 0;
}


// power_expansion 类型函数及其导数
int Non_Gaussianity_power_expansion_n(arb_t res, const arb_t x, const slong order, slong prec)
{
    arb_t s,t,A,B,C,D,E;
    
    arb_init(s);
    arb_init(t);
    arb_init(A);
    arb_init(B);
    arb_init(C);
    arb_init(D);
    arb_init(E);
    
    //power_expansion 非高斯性的函数为 ζ(r)=ζ_G(r) + 3/5 * f_NL * ζ_G^2 + 9/25*g_NL*ζ_G^3 + C*ζ_G^4 + D*ζ_G^5 + E*ζ_G^6 + ···
    //设常数 A=3/5*f_NL  常数 B=9/25*g_NL
    //y=x+A*x^2+B*x^3+C*x^4+D*x^5+E*x^6+F*x^7+G*x^8+H*x^9+I*x^10+J*x^11+K*x^12+L*x^13+M*x^14+N*x^15+O*x^16+P*x^17+Q*x^18+R*x^19+S*x^20+T*x^21
    
    arb_one(A);
    arb_mul_ui(A,A,3,prec);
    arb_div_ui(A,A,5,prec);
    arb_mul(A,A,Power_expansion_f,prec);
    
    arb_one(B);
    arb_mul_ui(B,B,9,prec);
    arb_div_ui(B,B,25,prec);
    arb_mul(B,B,Power_expansion_g,prec);
    
    arb_set(C,Power_expansion_four);
    arb_set(D,Power_expansion_five);
    arb_set(E,Power_expansion_six);
    
    switch(order)
    {
        case 0: //零阶导
            //y=x+A*x^2+B*x^3+C*x^4+D*x^5+E*x^6+F*x^7+G*x^8+H*x^9+I*x^10+J*x^11+K*x^12+L*x^13+M*x^14+N*x^15+O*x^16+P*x^17+Q*x^18+R*x^19+S*x^20+T*x^21
            
            arb_sqr(t,x,prec); //x+A*x^2
            arb_mul(t,t,A,prec);
            arb_add(t,t,x,prec);
            
            arb_pow_ui(s,x,3,prec);//B*x^3
            arb_mul(s,s,B,prec);
            arb_add(t,t,s,prec);
            
            arb_pow_ui(s,x,4,prec);//C*x^4
            arb_mul(s,s,C,prec);
            arb_add(t,t,s,prec);
            
            arb_pow_ui(s,x,5,prec);//D*x^5
            arb_mul(s,s,D,prec);
            arb_add(t,t,s,prec);
            
            arb_pow_ui(s,x,6,prec);//E*x^6
            arb_mul(s,s,E,prec);
            
            arb_add(res,t,s,prec);
            
            break;
        case 1: //一阶导
            //y=1+2*A*x+3*B*x^2+4*C*x^3+5*D*x^4+6*E*x^5+7*F*x^6+8*G*x^7+9*H*x^8+10*I*x^9+11*J*x^10+12*K*x^11+13*L*x^12+14*M*x^13+15*N*x^14+16*O*x^15+17*P*x^16+18*Q*x^17+19*R*x^18+20*S*x^19+21*T*x^20
            arb_mul_ui(t,x,2,prec);//1+A*2*x
            arb_mul(t,t,A,prec);
            arb_add_ui(t,t,1,prec);
            
            arb_sqr(s,x,prec); //B*3*x^2
            arb_mul_ui(s,s,3,prec);
            arb_mul(s,s,B,prec);
            arb_add(t,t,s,prec);
            
            arb_pow_ui(s,x,3,prec); //C*4*x^3
            arb_mul_ui(s,s,4,prec);
            arb_mul(s,s,C,prec);
            arb_add(t,t,s,prec);
            
            arb_pow_ui(s,x,4,prec); //D*5*x^4
            arb_mul_ui(s,s,5,prec);
            arb_mul(s,s,D,prec);
            arb_add(t,t,s,prec);
            
            arb_pow_ui(s,x,5,prec); //E*6*x^5
            arb_mul_ui(s,s,6,prec);
            arb_mul(s,s,E,prec);
            
            arb_add(res,t,s,prec);
            
            break;
        case 2: //二阶导
            //y=2*A+6*B*x+12*C*x^2+20*D*x^3+30*E*x^4+42*F*x^5+56*G*x^6+72*H*x^7+90*I*x^8+110*J*x^9+132*K*x^10+156*L*x^11+182*M*x^12+210*N*x^13+240*O*x^14+272*P*x^15+306*Q*x^16+342*R*x^17+380*S*x^18+420*T*x^19
            arb_mul_ui(t,A,2,prec); //2*A
            
            arb_mul_ui(s,x,6,prec); //B*6*x
            arb_mul(s,s,B,prec);
            arb_add(t,t,s,prec);
            
            arb_sqr(s,x,prec); //C*12*x^2
            arb_mul_ui(s,s,12,prec);
            arb_mul(s,s,C,prec);
            arb_add(t,t,s,prec);
            
            arb_pow_ui(s,x,3,prec); //D*20*x^3
            arb_mul_ui(s,s,20,prec);
            arb_mul(s,s,D,prec);
            arb_add(t,t,s,prec);
            
            arb_pow_ui(s,x,4,prec); //E*30*x^4
            arb_mul_ui(s,s,30,prec);
            arb_mul(s,s,E,prec);
            
            arb_add(res,t,s,prec);
            
            break;
        case 3: //三阶导
            //y=6*B+24*C*x+60*D*x^2+120*E*x^3+210*F*x^4+336*G*x^5+504*H*x^6+720*I*x^7+990*J*x^8+1320*K*x^9+1716*L*x^10+2184*M*x^11+2730*N*x^12+3360*O*x^13+4080*P*x^14+4896*Q*x^15+5814*R*x^16+6840*S*x^17+7980*T*x^18
            arb_mul_ui(t,B,6,prec); //B*6
            
            arb_mul_ui(s,x,24,prec); //C*24*x
            arb_mul(s,s,C,prec);
            arb_add(t,t,s,prec);
            
            arb_sqr(s,x,prec); //D*60*x^2
            arb_mul_ui(s,s,60,prec);
            arb_mul(s,s,D,prec);
            arb_add(t,t,s,prec);
            
            arb_pow_ui(s,x,3,prec); //E*120*x^3
            arb_mul_ui(s,s,120,prec);
            arb_mul(s,s,E,prec);
            
            arb_add(res,t,s,prec);
            
            break;
        case 4: //四阶导
            //y=24*C+120*D*x+360*E*x^2+840*F*x^3+1680*G*x^4+3024*H*x^5+5040*I*x^6+7920*J*x^7+11880*K*x^8+17160*L*x^9+24024*M*x^10+32760*N*x^11+43680*O*x^12+57120*P*x^13+73440*Q*x^14+93024*R*x^15+116280*S*x^16+143640*T*x^17
            arb_mul_ui(t,C,24,prec); //C*24
            
            arb_mul_ui(s,x,120,prec); //D*120*x
            arb_mul(s,s,D,prec);
            arb_add(t,t,s,prec);
            
            arb_sqr(s,x,prec); //E*360*x^2
            arb_mul_ui(s,s,360,prec);
            arb_mul(s,s,E,prec);
            
            arb_add(res,t,s,prec);
            
            break;
        default:
            printf("General -> non_gaussianity -> Non_Gaussianity_power_expansion_n 阶数n输入有误\n");
            exit(1);
    }
    
    
    arb_clear(s);
    arb_clear(t);
    arb_clear(A);
    arb_clear(B);
    arb_clear(C);
    arb_clear(D);
    arb_clear(E);
    
    return 0;
}



