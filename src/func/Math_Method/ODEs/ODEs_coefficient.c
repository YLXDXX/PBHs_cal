#include "../ODEs.h" 
#include <stdlib.h> 

//求于求解常微分方程的各种系数
//这里，作为全局变量存在

//用于误差估算所需系数
arb_t ODEs_step_max_factor,ODEs_step_min_factor,ODEs_step_factor;
arb_ptr ODEs_step_factor_help;//此变量无实际作用，仅用于判断初始化是否完成

void ODEs_get_step_cal_coe(slong prec)
{
    if(ODEs_step_factor_help != NULL) //判断是否已设定系数
    {
        return;
    }
    
    arb_init(ODEs_step_max_factor); //用于自动步长
    arb_init(ODEs_step_min_factor);
    arb_init(ODEs_step_factor);
    
    arb_set_str(ODEs_step_max_factor,"5",prec);
    arb_set_str(ODEs_step_min_factor,"0.2",prec);
    arb_set_str(ODEs_step_factor,"0.9",prec);
    
    ODEs_step_factor_help=_arb_vec_init(1); //用于判断初始化是否完成
    arb_one(ODEs_step_factor_help);
}


// RFK45 方法
//求k_1 .. k_6 相关系数
arb_ptr ODEs_RFK45_coe_ah;
arb_ptr ODEs_RFK45_coe_b1,ODEs_RFK45_coe_b2;
arb_ptr ODEs_RFK45_coe_b3,ODEs_RFK45_coe_b4,ODEs_RFK45_coe_b5;
arb_ptr ODEs_RFK45_coe_ck_s,ODEs_RFK45_coe_ck_t;

void ODEs_get_RFK45_cal_coe(slong prec)
{
    if(ODEs_RFK45_coe_ah != NULL) //判断是否已设定系数
    {
        return;
    }
    
    ODEs_RFK45_coe_ah=_arb_vec_init(7);
    ODEs_RFK45_coe_b1=_arb_vec_init(7);
    ODEs_RFK45_coe_b2=_arb_vec_init(7);
    ODEs_RFK45_coe_b3=_arb_vec_init(7);
    ODEs_RFK45_coe_b4=_arb_vec_init(7);
    ODEs_RFK45_coe_b5=_arb_vec_init(7);
    
    ODEs_RFK45_coe_ck_s=_arb_vec_init(7);
    ODEs_RFK45_coe_ck_t=_arb_vec_init(7);
    
    for(int i=0; i<7; i++ )
    {
        arb_one(ODEs_RFK45_coe_ah+i);
        arb_one(ODEs_RFK45_coe_b1+i);
        arb_one(ODEs_RFK45_coe_b2+i);
        arb_one(ODEs_RFK45_coe_b3+i);
        arb_one(ODEs_RFK45_coe_b4+i);
        arb_one(ODEs_RFK45_coe_b5+i);
        arb_one(ODEs_RFK45_coe_ck_s+i);
        arb_one(ODEs_RFK45_coe_ck_t+i);
    }
    
    //ODEs_RFK45_coe_ah k_2 ... k_6
    arb_zero(ODEs_RFK45_coe_ah+0);
    arb_zero(ODEs_RFK45_coe_ah+1);
    
    arb_div_ui(ODEs_RFK45_coe_ah+2,ODEs_RFK45_coe_ah+2,4,prec); //1/4
    
    arb_mul_ui(ODEs_RFK45_coe_ah+3,ODEs_RFK45_coe_ah+3,3,prec); //3/8
    arb_div_ui(ODEs_RFK45_coe_ah+3,ODEs_RFK45_coe_ah+3,8,prec);
    
    arb_mul_ui(ODEs_RFK45_coe_ah+4,ODEs_RFK45_coe_ah+4,12,prec); //12/13
    arb_div_ui(ODEs_RFK45_coe_ah+4,ODEs_RFK45_coe_ah+4,13,prec);
    
    arb_one(ODEs_RFK45_coe_ah+5); //1
    
    arb_div_ui(ODEs_RFK45_coe_ah+6,ODEs_RFK45_coe_ah+6,2,prec); //1/2
    
    
    //ODEs_RFK45_coe_b1 k_2 ... k_6
    arb_zero(ODEs_RFK45_coe_b1+0);
    arb_zero(ODEs_RFK45_coe_b1+1);
    
    arb_div_ui(ODEs_RFK45_coe_b1+2,ODEs_RFK45_coe_b1+2,4,prec); //1/4
    
    arb_mul_ui(ODEs_RFK45_coe_b1+3,ODEs_RFK45_coe_b1+3,3,prec); //3/32
    arb_div_ui(ODEs_RFK45_coe_b1+3,ODEs_RFK45_coe_b1+3,32,prec);
    
    arb_mul_ui(ODEs_RFK45_coe_b1+4,ODEs_RFK45_coe_b1+4,1932,prec); //1932/2197
    arb_div_ui(ODEs_RFK45_coe_b1+4,ODEs_RFK45_coe_b1+4,2197,prec);
    
    arb_mul_ui(ODEs_RFK45_coe_b1+5,ODEs_RFK45_coe_b1+5,439,prec); //439/216
    arb_div_ui(ODEs_RFK45_coe_b1+5,ODEs_RFK45_coe_b1+5,216,prec);
    
    arb_mul_si(ODEs_RFK45_coe_b1+6,ODEs_RFK45_coe_b1+6,-8,prec); //-8/27
    arb_div_ui(ODEs_RFK45_coe_b1+6,ODEs_RFK45_coe_b1+6,27,prec);
    
    
    //ODEs_RFK45_coe_b2 k_3 ... k_6
    arb_zero(ODEs_RFK45_coe_b2+0);
    arb_zero(ODEs_RFK45_coe_b2+1);
    arb_zero(ODEs_RFK45_coe_b2+2);
    
    arb_mul_ui(ODEs_RFK45_coe_b2+3,ODEs_RFK45_coe_b2+3,9,prec); //9/32
    arb_div_ui(ODEs_RFK45_coe_b2+3,ODEs_RFK45_coe_b2+3,32,prec);
    
    arb_mul_si(ODEs_RFK45_coe_b2+4,ODEs_RFK45_coe_b2+4,-7200,prec); //-7200/2197
    arb_div_ui(ODEs_RFK45_coe_b2+4,ODEs_RFK45_coe_b2+4,2197,prec);
    
    arb_mul_si(ODEs_RFK45_coe_b2+5,ODEs_RFK45_coe_b2+5,-8,prec); //-8
    
    arb_mul_si(ODEs_RFK45_coe_b2+6,ODEs_RFK45_coe_b2+6,2,prec); //2
    
    
    //ODEs_RFK45_coe_b3 k_4 ... k_6
    arb_zero(ODEs_RFK45_coe_b3+0);
    arb_zero(ODEs_RFK45_coe_b3+1);
    arb_zero(ODEs_RFK45_coe_b3+2);
    arb_zero(ODEs_RFK45_coe_b3+3);
    
    arb_mul_ui(ODEs_RFK45_coe_b3+4,ODEs_RFK45_coe_b3+4,7296,prec); //7296/2197
    arb_div_ui(ODEs_RFK45_coe_b3+4,ODEs_RFK45_coe_b3+4,2197,prec);
    
    arb_mul_si(ODEs_RFK45_coe_b3+5,ODEs_RFK45_coe_b3+5,3680,prec); //3680/513
    arb_div_ui(ODEs_RFK45_coe_b3+5,ODEs_RFK45_coe_b3+5,513,prec);
    
    arb_mul_si(ODEs_RFK45_coe_b3+6,ODEs_RFK45_coe_b3+6,-3544,prec); //-3544/2565
    arb_div_ui(ODEs_RFK45_coe_b3+6,ODEs_RFK45_coe_b3+6,2565,prec);
    
    //ODEs_RFK45_coe_b4 k_5 ... k_6
    arb_zero(ODEs_RFK45_coe_b4+0);
    arb_zero(ODEs_RFK45_coe_b4+1);
    arb_zero(ODEs_RFK45_coe_b4+2);
    arb_zero(ODEs_RFK45_coe_b4+3);
    arb_zero(ODEs_RFK45_coe_b4+4);
    
    arb_mul_si(ODEs_RFK45_coe_b4+5,ODEs_RFK45_coe_b4+5,-845,prec); //-845/4104
    arb_div_ui(ODEs_RFK45_coe_b4+5,ODEs_RFK45_coe_b4+5,4104,prec);
    
    arb_mul_si(ODEs_RFK45_coe_b4+6,ODEs_RFK45_coe_b4+6,1859,prec); //1859/4104
    arb_div_ui(ODEs_RFK45_coe_b4+6,ODEs_RFK45_coe_b4+6,4104,prec);
    
    
    //ODEs_RFK45_coe_b5 k_6 ... k_6
    arb_zero(ODEs_RFK45_coe_b5+0);
    arb_zero(ODEs_RFK45_coe_b5+1);
    arb_zero(ODEs_RFK45_coe_b5+2);
    arb_zero(ODEs_RFK45_coe_b5+3);
    arb_zero(ODEs_RFK45_coe_b5+4);
    arb_zero(ODEs_RFK45_coe_b5+5);
    
    arb_mul_si(ODEs_RFK45_coe_b5+6,ODEs_RFK45_coe_b5+6,-11,prec); //-11/40
    arb_div_ui(ODEs_RFK45_coe_b5+6,ODEs_RFK45_coe_b5+6,40,prec);
    
    
    //ODEs_RFK45_coe_ck_s  k_1 ... k_5
    arb_zero(ODEs_RFK45_coe_ck_s+0);
    
    arb_mul_si(ODEs_RFK45_coe_ck_s+1,ODEs_RFK45_coe_ck_s+1,25,prec); //25/216
    arb_div_ui(ODEs_RFK45_coe_ck_s+1,ODEs_RFK45_coe_ck_s+1,216,prec);
    
    arb_zero(ODEs_RFK45_coe_ck_s+2); //0
    
    arb_mul_si(ODEs_RFK45_coe_ck_s+3,ODEs_RFK45_coe_ck_s+3,1408,prec); //1408/2565
    arb_div_ui(ODEs_RFK45_coe_ck_s+3,ODEs_RFK45_coe_ck_s+3,2565,prec);
    
    arb_mul_si(ODEs_RFK45_coe_ck_s+4,ODEs_RFK45_coe_ck_s+4,2197,prec); //2197/4104
    arb_div_ui(ODEs_RFK45_coe_ck_s+4,ODEs_RFK45_coe_ck_s+4,4104,prec);
    
    arb_mul_si(ODEs_RFK45_coe_ck_s+5,ODEs_RFK45_coe_ck_s+5,-1,prec); //-1/5
    arb_div_ui(ODEs_RFK45_coe_ck_s+5,ODEs_RFK45_coe_ck_s+5,5,prec);
    
    arb_zero(ODEs_RFK45_coe_ck_s+6); //0
    
    
    //ODEs_RFK45_coe_ck_t  k_1 ... k_5
    arb_zero(ODEs_RFK45_coe_ck_t+0);
    
    arb_mul_si(ODEs_RFK45_coe_ck_t+1,ODEs_RFK45_coe_ck_t+1,16,prec); //16/135
    arb_div_ui(ODEs_RFK45_coe_ck_t+1,ODEs_RFK45_coe_ck_t+1,135,prec);
    
    arb_zero(ODEs_RFK45_coe_ck_t+2); //0
    
    arb_mul_si(ODEs_RFK45_coe_ck_t+3,ODEs_RFK45_coe_ck_t+3,6656,prec); //6656/12825
    arb_div_ui(ODEs_RFK45_coe_ck_t+3,ODEs_RFK45_coe_ck_t+3,12825,prec);
    
    arb_mul_si(ODEs_RFK45_coe_ck_t+4,ODEs_RFK45_coe_ck_t+4,28561,prec); //28561/56430
    arb_div_ui(ODEs_RFK45_coe_ck_t+4,ODEs_RFK45_coe_ck_t+4,56430,prec);
    
    arb_mul_si(ODEs_RFK45_coe_ck_t+5,ODEs_RFK45_coe_ck_t+5,-9,prec); //-9/50
    arb_div_ui(ODEs_RFK45_coe_ck_t+5,ODEs_RFK45_coe_ck_t+5,50,prec);
    
    arb_mul_si(ODEs_RFK45_coe_ck_t+6,ODEs_RFK45_coe_ck_t+6,2,prec); //2/55
    arb_div_ui(ODEs_RFK45_coe_ck_t+6,ODEs_RFK45_coe_ck_t+6,55,prec);
    
}



// DOPRI54 方法，所需各系数
//其中 dense output 参见 https://math.stackexchange.com/questions/2947231/how-can-i-derive-the-dense-output-of-ode45
arb_ptr ODEs_DOPRI54_coe_c;
arb_ptr ODEs_DOPRI54_coe_b_up;
arb_ptr ODEs_DOPRI54_coe_b_dw;
arb_ptr ODEs_DOPRI54_coe_a_1;
arb_ptr ODEs_DOPRI54_coe_a_2;
arb_ptr ODEs_DOPRI54_coe_a_3;
arb_ptr ODEs_DOPRI54_coe_a_4;
arb_ptr ODEs_DOPRI54_coe_a_5;
arb_ptr ODEs_DOPRI54_coe_a_6;
arb_ptr ODEs_DOPRI54_coe_dense_out;

void ODEs_get_DOPRI54_cal_coe(slong prec)
{
    if(ODEs_DOPRI54_coe_c != NULL) //判断是否已设定系数
    {
        return;
    }
    
    ODEs_DOPRI54_coe_c=_arb_vec_init(7);
    ODEs_DOPRI54_coe_b_up=_arb_vec_init(7);
    ODEs_DOPRI54_coe_b_dw=_arb_vec_init(7);
    ODEs_DOPRI54_coe_a_1=_arb_vec_init(7);
    ODEs_DOPRI54_coe_a_2=_arb_vec_init(7);
    ODEs_DOPRI54_coe_a_3=_arb_vec_init(7);
    ODEs_DOPRI54_coe_a_4=_arb_vec_init(7);
    ODEs_DOPRI54_coe_a_5=_arb_vec_init(7);
    ODEs_DOPRI54_coe_a_6=_arb_vec_init(7);
    ODEs_DOPRI54_coe_dense_out=_arb_vec_init(7);
    
    //所有元素设为零
    _arb_vec_zero(ODEs_DOPRI54_coe_c, 7);
    _arb_vec_zero(ODEs_DOPRI54_coe_b_up, 7);
    _arb_vec_zero(ODEs_DOPRI54_coe_b_dw, 7);
    _arb_vec_zero(ODEs_DOPRI54_coe_a_1, 7);
    _arb_vec_zero(ODEs_DOPRI54_coe_a_2, 7);
    _arb_vec_zero(ODEs_DOPRI54_coe_a_3, 7);
    _arb_vec_zero(ODEs_DOPRI54_coe_a_4, 7);
    _arb_vec_zero(ODEs_DOPRI54_coe_a_5, 7);
    _arb_vec_zero(ODEs_DOPRI54_coe_a_6, 7);
    _arb_vec_zero(ODEs_DOPRI54_coe_dense_out, 7);
    
    //赋值
    arb_set_si(ODEs_DOPRI54_coe_c+1,1);
    arb_div_si(ODEs_DOPRI54_coe_c+1,ODEs_DOPRI54_coe_c+1,5,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_c+2,3);
    arb_div_si(ODEs_DOPRI54_coe_c+2,ODEs_DOPRI54_coe_c+2,10,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_c+3,4);
    arb_div_si(ODEs_DOPRI54_coe_c+3,ODEs_DOPRI54_coe_c+3,5,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_c+4,8);
    arb_div_si(ODEs_DOPRI54_coe_c+4,ODEs_DOPRI54_coe_c+4,9,prec);
    
    arb_one(ODEs_DOPRI54_coe_c+5);
    arb_one(ODEs_DOPRI54_coe_c+6);
    
    
    arb_set_si(ODEs_DOPRI54_coe_b_up+0,35);
    arb_div_si(ODEs_DOPRI54_coe_b_up+0,ODEs_DOPRI54_coe_b_up+0,384,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_b_up+2,500);
    arb_div_si(ODEs_DOPRI54_coe_b_up+2,ODEs_DOPRI54_coe_b_up+2,1113,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_b_up+3,125);
    arb_div_si(ODEs_DOPRI54_coe_b_up+3,ODEs_DOPRI54_coe_b_up+3,192,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_b_up+4,-2187);
    arb_div_si(ODEs_DOPRI54_coe_b_up+4,ODEs_DOPRI54_coe_b_up+4,6784,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_b_up+5,11);
    arb_div_si(ODEs_DOPRI54_coe_b_up+5,ODEs_DOPRI54_coe_b_up+5,84,prec);
    
    
    arb_set_si(ODEs_DOPRI54_coe_b_dw+0,5179);
    arb_div_si(ODEs_DOPRI54_coe_b_dw+0,ODEs_DOPRI54_coe_b_dw+0,57600,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_b_dw+2,7571);
    arb_div_si(ODEs_DOPRI54_coe_b_dw+2,ODEs_DOPRI54_coe_b_dw+2,16695,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_b_dw+3,393);
    arb_div_si(ODEs_DOPRI54_coe_b_dw+3,ODEs_DOPRI54_coe_b_dw+3,640,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_b_dw+4,-92097);
    arb_div_si(ODEs_DOPRI54_coe_b_dw+4,ODEs_DOPRI54_coe_b_dw+4,339200,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_b_dw+5,187);
    arb_div_si(ODEs_DOPRI54_coe_b_dw+5,ODEs_DOPRI54_coe_b_dw+5,2100,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_b_dw+6,1);
    arb_div_si(ODEs_DOPRI54_coe_b_dw+6,ODEs_DOPRI54_coe_b_dw+6,40,prec);
    
    
    
    arb_set_si(ODEs_DOPRI54_coe_a_1+0,1);
    arb_div_si(ODEs_DOPRI54_coe_a_1+0,ODEs_DOPRI54_coe_a_1+0,5,prec);
    
    
    arb_set_si(ODEs_DOPRI54_coe_a_2+0,3);
    arb_div_si(ODEs_DOPRI54_coe_a_2+0,ODEs_DOPRI54_coe_a_2+0,40,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_a_2+1,9);
    arb_div_si(ODEs_DOPRI54_coe_a_2+1,ODEs_DOPRI54_coe_a_2+1,40,prec);
    
    
    arb_set_si(ODEs_DOPRI54_coe_a_3+0,44);
    arb_div_si(ODEs_DOPRI54_coe_a_3+0,ODEs_DOPRI54_coe_a_3+0,45,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_a_3+1,-56);
    arb_div_si(ODEs_DOPRI54_coe_a_3+1,ODEs_DOPRI54_coe_a_3+1,15,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_a_3+2,32);
    arb_div_si(ODEs_DOPRI54_coe_a_3+2,ODEs_DOPRI54_coe_a_3+2,9,prec);
    
    
    arb_set_si(ODEs_DOPRI54_coe_a_4+0,19372);
    arb_div_si(ODEs_DOPRI54_coe_a_4+0,ODEs_DOPRI54_coe_a_4+0,6561,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_a_4+1,-25360);
    arb_div_si(ODEs_DOPRI54_coe_a_4+1,ODEs_DOPRI54_coe_a_4+1,2187,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_a_4+2,64448);
    arb_div_si(ODEs_DOPRI54_coe_a_4+2,ODEs_DOPRI54_coe_a_4+2,6561,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_a_4+3,-212);
    arb_div_si(ODEs_DOPRI54_coe_a_4+3,ODEs_DOPRI54_coe_a_4+3,729,prec);
    
    
    arb_set_si(ODEs_DOPRI54_coe_a_5+0,9017);
    arb_div_si(ODEs_DOPRI54_coe_a_5+0,ODEs_DOPRI54_coe_a_5+0,3168,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_a_5+1,-355);
    arb_div_si(ODEs_DOPRI54_coe_a_5+1,ODEs_DOPRI54_coe_a_5+1,33,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_a_5+2,46732);
    arb_div_si(ODEs_DOPRI54_coe_a_5+2,ODEs_DOPRI54_coe_a_5+2,5247,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_a_5+3,49);
    arb_div_si(ODEs_DOPRI54_coe_a_5+3,ODEs_DOPRI54_coe_a_5+3,176,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_a_5+4,-5103);
    arb_div_si(ODEs_DOPRI54_coe_a_5+4,ODEs_DOPRI54_coe_a_5+4,18656,prec);
    
    
    arb_set_si(ODEs_DOPRI54_coe_a_6+0,35);
    arb_div_si(ODEs_DOPRI54_coe_a_6+0,ODEs_DOPRI54_coe_a_6+0,384,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_a_6+2,500);
    arb_div_si(ODEs_DOPRI54_coe_a_6+2,ODEs_DOPRI54_coe_a_6+2,1113,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_a_6+3,125);
    arb_div_si(ODEs_DOPRI54_coe_a_6+3,ODEs_DOPRI54_coe_a_6+3,192,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_a_6+4,-2187);
    arb_div_si(ODEs_DOPRI54_coe_a_6+4,ODEs_DOPRI54_coe_a_6+4,6784,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_a_6+5,11);
    arb_div_si(ODEs_DOPRI54_coe_a_6+5,ODEs_DOPRI54_coe_a_6+5,84,prec);
    
    
    arb_set_si(ODEs_DOPRI54_coe_dense_out+0,-12715105075);
    arb_div_si(ODEs_DOPRI54_coe_dense_out+0,ODEs_DOPRI54_coe_dense_out+0,11282082432,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_dense_out+2,87487479700);
    arb_div_si(ODEs_DOPRI54_coe_dense_out+2,ODEs_DOPRI54_coe_dense_out+2,32700410799,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_dense_out+3,-10690763975);
    arb_div_si(ODEs_DOPRI54_coe_dense_out+3,ODEs_DOPRI54_coe_dense_out+3,1880347072,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_dense_out+4,701980252875);
    arb_div_si(ODEs_DOPRI54_coe_dense_out+4,ODEs_DOPRI54_coe_dense_out+4,199316789632,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_dense_out+5,-1453857185);
    arb_div_si(ODEs_DOPRI54_coe_dense_out+5,ODEs_DOPRI54_coe_dense_out+5,822651844,prec);
    
    arb_set_si(ODEs_DOPRI54_coe_dense_out+6,69997945);
    arb_div_si(ODEs_DOPRI54_coe_dense_out+6,ODEs_DOPRI54_coe_dense_out+6,29380423,prec);
    
}

// DOP853 方法，所需各系数
// 该系数表来源于 https://github.com/SciML/OrdinaryDiffEq.jl/blob/master/lib/OrdinaryDiffEqHighOrderRK/src/high_order_rk_tableaus.jl
arb_t ODEs_DOP853_coe_c7, ODEs_DOP853_coe_c8, ODEs_DOP853_coe_c9, ODEs_DOP853_coe_c10, ODEs_DOP853_coe_c11, ODEs_DOP853_coe_c6, ODEs_DOP853_coe_c5, ODEs_DOP853_coe_c4, ODEs_DOP853_coe_c3, ODEs_DOP853_coe_c2, ODEs_DOP853_coe_b1, ODEs_DOP853_coe_b6, ODEs_DOP853_coe_b7, ODEs_DOP853_coe_b8, ODEs_DOP853_coe_b9, ODEs_DOP853_coe_b10, ODEs_DOP853_coe_b11;
arb_t ODEs_DOP853_coe_b12, ODEs_DOP853_coe_btilde1, ODEs_DOP853_coe_btilde6, ODEs_DOP853_coe_btilde7, ODEs_DOP853_coe_btilde8, ODEs_DOP853_coe_btilde9, ODEs_DOP853_coe_btilde10, ODEs_DOP853_coe_btilde11;
arb_t ODEs_DOP853_coe_btilde12, ODEs_DOP853_coe_er1, ODEs_DOP853_coe_er6, ODEs_DOP853_coe_er7, ODEs_DOP853_coe_er8, ODEs_DOP853_coe_er9, ODEs_DOP853_coe_er10, ODEs_DOP853_coe_er11, ODEs_DOP853_coe_er12, ODEs_DOP853_coe_a0201, ODEs_DOP853_coe_a0301;
arb_t ODEs_DOP853_coe_a0302, ODEs_DOP853_coe_a0401, ODEs_DOP853_coe_a0403, ODEs_DOP853_coe_a0501, ODEs_DOP853_coe_a0503, ODEs_DOP853_coe_a0504, ODEs_DOP853_coe_a0601, ODEs_DOP853_coe_a0604, ODEs_DOP853_coe_a0605, ODEs_DOP853_coe_a0701;
arb_t ODEs_DOP853_coe_a0704, ODEs_DOP853_coe_a0705, ODEs_DOP853_coe_a0706, ODEs_DOP853_coe_a0801, ODEs_DOP853_coe_a0804, ODEs_DOP853_coe_a0805, ODEs_DOP853_coe_a0806, ODEs_DOP853_coe_a0807, ODEs_DOP853_coe_a0901, ODEs_DOP853_coe_a0904;
arb_t ODEs_DOP853_coe_a0905, ODEs_DOP853_coe_a0906, ODEs_DOP853_coe_a0907, ODEs_DOP853_coe_a0908, ODEs_DOP853_coe_a1001, ODEs_DOP853_coe_a1004, ODEs_DOP853_coe_a1005, ODEs_DOP853_coe_a1006, ODEs_DOP853_coe_a1007, ODEs_DOP853_coe_a1008;
arb_t ODEs_DOP853_coe_a1009, ODEs_DOP853_coe_a1101, ODEs_DOP853_coe_a1104, ODEs_DOP853_coe_a1105, ODEs_DOP853_coe_a1106, ODEs_DOP853_coe_a1107, ODEs_DOP853_coe_a1108, ODEs_DOP853_coe_a1109, ODEs_DOP853_coe_a1110, ODEs_DOP853_coe_a1201;
arb_t ODEs_DOP853_coe_a1204, ODEs_DOP853_coe_a1205, ODEs_DOP853_coe_a1206, ODEs_DOP853_coe_a1207, ODEs_DOP853_coe_a1208, ODEs_DOP853_coe_a1209, ODEs_DOP853_coe_a1210, ODEs_DOP853_coe_a1211, ODEs_DOP853_coe_c14, ODEs_DOP853_coe_c15, ODEs_DOP853_coe_c16;
arb_t ODEs_DOP853_coe_a1401, ODEs_DOP853_coe_a1407, ODEs_DOP853_coe_a1408, ODEs_DOP853_coe_a1409, ODEs_DOP853_coe_a1410, ODEs_DOP853_coe_a1411, ODEs_DOP853_coe_a1412, ODEs_DOP853_coe_a1413, ODEs_DOP853_coe_a1501, ODEs_DOP853_coe_a1506;
arb_t ODEs_DOP853_coe_a1507, ODEs_DOP853_coe_a1508, ODEs_DOP853_coe_a1511, ODEs_DOP853_coe_a1512, ODEs_DOP853_coe_a1513, ODEs_DOP853_coe_a1514, ODEs_DOP853_coe_a1601, ODEs_DOP853_coe_a1606, ODEs_DOP853_coe_a1607, ODEs_DOP853_coe_a1608;
arb_t ODEs_DOP853_coe_a1609, ODEs_DOP853_coe_a1613, ODEs_DOP853_coe_a1614, ODEs_DOP853_coe_a1615, ODEs_DOP853_coe_d401, ODEs_DOP853_coe_d406, ODEs_DOP853_coe_d407, ODEs_DOP853_coe_d408, ODEs_DOP853_coe_d409, ODEs_DOP853_coe_d410, ODEs_DOP853_coe_d411;
arb_t ODEs_DOP853_coe_d412, ODEs_DOP853_coe_d413, ODEs_DOP853_coe_d414, ODEs_DOP853_coe_d415, ODEs_DOP853_coe_d416, ODEs_DOP853_coe_d501, ODEs_DOP853_coe_d506, ODEs_DOP853_coe_d507, ODEs_DOP853_coe_d508, ODEs_DOP853_coe_d509, ODEs_DOP853_coe_d510, ODEs_DOP853_coe_d511;
arb_t ODEs_DOP853_coe_d512, ODEs_DOP853_coe_d513, ODEs_DOP853_coe_d514, ODEs_DOP853_coe_d515, ODEs_DOP853_coe_d516, ODEs_DOP853_coe_d601, ODEs_DOP853_coe_d606, ODEs_DOP853_coe_d607, ODEs_DOP853_coe_d608, ODEs_DOP853_coe_d609, ODEs_DOP853_coe_d610, ODEs_DOP853_coe_d611;
arb_t ODEs_DOP853_coe_d612, ODEs_DOP853_coe_d613, ODEs_DOP853_coe_d614, ODEs_DOP853_coe_d615, ODEs_DOP853_coe_d616, ODEs_DOP853_coe_d701, ODEs_DOP853_coe_d706, ODEs_DOP853_coe_d707, ODEs_DOP853_coe_d708, ODEs_DOP853_coe_d709, ODEs_DOP853_coe_d710, ODEs_DOP853_coe_d711;
arb_t ODEs_DOP853_coe_d712, ODEs_DOP853_coe_d713, ODEs_DOP853_coe_d714, ODEs_DOP853_coe_d715, ODEs_DOP853_coe_d716;
arb_t  ODEs_DOP853_coe_bhh1,ODEs_DOP853_coe_bhh2,ODEs_DOP853_coe_bhh3;

arb_ptr ODEs_get_DOP853_cal_coe_help;//此变量无实际作用，仅用于判断初始化是否完成

void ODEs_get_DOP853_cal_coe(slong prec)
{
    if(ODEs_get_DOP853_cal_coe_help != NULL) //判断是否已设定系数
    {
        return;
    }
    
    arb_init(ODEs_DOP853_coe_c7);
    arb_init(ODEs_DOP853_coe_c8);
    arb_init(ODEs_DOP853_coe_c9);
    arb_init(ODEs_DOP853_coe_c10);
    arb_init(ODEs_DOP853_coe_c11);
    arb_init(ODEs_DOP853_coe_c6);
    arb_init(ODEs_DOP853_coe_c5);
    arb_init(ODEs_DOP853_coe_c4);
    arb_init(ODEs_DOP853_coe_c3);
    arb_init(ODEs_DOP853_coe_c2);
    arb_init(ODEs_DOP853_coe_b1);
    arb_init(ODEs_DOP853_coe_b6);
    arb_init(ODEs_DOP853_coe_b7);
    arb_init(ODEs_DOP853_coe_b8);
    arb_init(ODEs_DOP853_coe_b9);
    arb_init(ODEs_DOP853_coe_b10);
    arb_init(ODEs_DOP853_coe_b11);
    arb_init(ODEs_DOP853_coe_b12);
    
    arb_init(ODEs_DOP853_coe_bhh1);
    arb_init(ODEs_DOP853_coe_bhh2);
    arb_init(ODEs_DOP853_coe_bhh3);
    
    arb_init(ODEs_DOP853_coe_btilde1);
    arb_init(ODEs_DOP853_coe_btilde6);
    arb_init(ODEs_DOP853_coe_btilde7);
    arb_init(ODEs_DOP853_coe_btilde8);
    arb_init(ODEs_DOP853_coe_btilde9);
    arb_init(ODEs_DOP853_coe_btilde10);
    arb_init(ODEs_DOP853_coe_btilde11);
    arb_init(ODEs_DOP853_coe_btilde12);
    arb_init(ODEs_DOP853_coe_er1);
    arb_init(ODEs_DOP853_coe_er6);
    arb_init(ODEs_DOP853_coe_er7);
    arb_init(ODEs_DOP853_coe_er8);
    arb_init(ODEs_DOP853_coe_er9);
    arb_init(ODEs_DOP853_coe_er10);
    arb_init(ODEs_DOP853_coe_er11);
    arb_init(ODEs_DOP853_coe_er12);
    arb_init(ODEs_DOP853_coe_a0201);
    arb_init(ODEs_DOP853_coe_a0301);
    arb_init(ODEs_DOP853_coe_a0302);
    arb_init(ODEs_DOP853_coe_a0401);
    arb_init(ODEs_DOP853_coe_a0403);
    arb_init(ODEs_DOP853_coe_a0501);
    arb_init(ODEs_DOP853_coe_a0503);
    arb_init(ODEs_DOP853_coe_a0504);
    arb_init(ODEs_DOP853_coe_a0601);
    arb_init(ODEs_DOP853_coe_a0604);
    arb_init(ODEs_DOP853_coe_a0605);
    arb_init(ODEs_DOP853_coe_a0701);
    arb_init(ODEs_DOP853_coe_a0704);
    arb_init(ODEs_DOP853_coe_a0705);
    arb_init(ODEs_DOP853_coe_a0706);
    arb_init(ODEs_DOP853_coe_a0801);
    arb_init(ODEs_DOP853_coe_a0804);
    arb_init(ODEs_DOP853_coe_a0805);
    arb_init(ODEs_DOP853_coe_a0806);
    arb_init(ODEs_DOP853_coe_a0807);
    arb_init(ODEs_DOP853_coe_a0901);
    arb_init(ODEs_DOP853_coe_a0904);
    arb_init(ODEs_DOP853_coe_a0905);
    arb_init(ODEs_DOP853_coe_a0906);
    arb_init(ODEs_DOP853_coe_a0907);
    arb_init(ODEs_DOP853_coe_a0908);
    arb_init(ODEs_DOP853_coe_a1001);
    arb_init(ODEs_DOP853_coe_a1004);
    arb_init(ODEs_DOP853_coe_a1005);
    arb_init(ODEs_DOP853_coe_a1006);
    arb_init(ODEs_DOP853_coe_a1007);
    arb_init(ODEs_DOP853_coe_a1008);
    arb_init(ODEs_DOP853_coe_a1009);
    arb_init(ODEs_DOP853_coe_a1101);
    arb_init(ODEs_DOP853_coe_a1104);
    arb_init(ODEs_DOP853_coe_a1105);
    arb_init(ODEs_DOP853_coe_a1106);
    arb_init(ODEs_DOP853_coe_a1107);
    arb_init(ODEs_DOP853_coe_a1108);
    arb_init(ODEs_DOP853_coe_a1109);
    arb_init(ODEs_DOP853_coe_a1110);
    arb_init(ODEs_DOP853_coe_a1201);
    arb_init(ODEs_DOP853_coe_a1204);
    arb_init(ODEs_DOP853_coe_a1205);
    arb_init(ODEs_DOP853_coe_a1206);
    arb_init(ODEs_DOP853_coe_a1207);
    arb_init(ODEs_DOP853_coe_a1208);
    arb_init(ODEs_DOP853_coe_a1209);
    arb_init(ODEs_DOP853_coe_a1210);
    arb_init(ODEs_DOP853_coe_a1211);
    arb_init(ODEs_DOP853_coe_c14);
    arb_init(ODEs_DOP853_coe_c15);
    arb_init(ODEs_DOP853_coe_c16);
    arb_init(ODEs_DOP853_coe_a1401);
    arb_init(ODEs_DOP853_coe_a1407);
    arb_init(ODEs_DOP853_coe_a1408);
    arb_init(ODEs_DOP853_coe_a1409);
    arb_init(ODEs_DOP853_coe_a1410);
    arb_init(ODEs_DOP853_coe_a1411);
    arb_init(ODEs_DOP853_coe_a1412);
    arb_init(ODEs_DOP853_coe_a1413);
    arb_init(ODEs_DOP853_coe_a1501);
    arb_init(ODEs_DOP853_coe_a1506);
    arb_init(ODEs_DOP853_coe_a1507);
    arb_init(ODEs_DOP853_coe_a1508);
    arb_init(ODEs_DOP853_coe_a1511);
    arb_init(ODEs_DOP853_coe_a1512);
    arb_init(ODEs_DOP853_coe_a1513);
    arb_init(ODEs_DOP853_coe_a1514);
    arb_init(ODEs_DOP853_coe_a1601);
    arb_init(ODEs_DOP853_coe_a1606);
    arb_init(ODEs_DOP853_coe_a1607);
    arb_init(ODEs_DOP853_coe_a1608);
    arb_init(ODEs_DOP853_coe_a1609);
    arb_init(ODEs_DOP853_coe_a1613);
    arb_init(ODEs_DOP853_coe_a1614);
    arb_init(ODEs_DOP853_coe_a1615);
    arb_init(ODEs_DOP853_coe_d401);
    arb_init(ODEs_DOP853_coe_d406);
    arb_init(ODEs_DOP853_coe_d407);
    arb_init(ODEs_DOP853_coe_d408);
    arb_init(ODEs_DOP853_coe_d409);
    arb_init(ODEs_DOP853_coe_d410);
    arb_init(ODEs_DOP853_coe_d411);
    arb_init(ODEs_DOP853_coe_d412);
    arb_init(ODEs_DOP853_coe_d413);
    arb_init(ODEs_DOP853_coe_d414);
    arb_init(ODEs_DOP853_coe_d415);
    arb_init(ODEs_DOP853_coe_d416);
    arb_init(ODEs_DOP853_coe_d501);
    arb_init(ODEs_DOP853_coe_d506);
    arb_init(ODEs_DOP853_coe_d507);
    arb_init(ODEs_DOP853_coe_d508);
    arb_init(ODEs_DOP853_coe_d509);
    arb_init(ODEs_DOP853_coe_d510);
    arb_init(ODEs_DOP853_coe_d511);
    arb_init(ODEs_DOP853_coe_d512);
    arb_init(ODEs_DOP853_coe_d513);
    arb_init(ODEs_DOP853_coe_d514);
    arb_init(ODEs_DOP853_coe_d515);
    arb_init(ODEs_DOP853_coe_d516);
    arb_init(ODEs_DOP853_coe_d601);
    arb_init(ODEs_DOP853_coe_d606);
    arb_init(ODEs_DOP853_coe_d607);
    arb_init(ODEs_DOP853_coe_d608);
    arb_init(ODEs_DOP853_coe_d609);
    arb_init(ODEs_DOP853_coe_d610);
    arb_init(ODEs_DOP853_coe_d611);
    arb_init(ODEs_DOP853_coe_d612);
    arb_init(ODEs_DOP853_coe_d613);
    arb_init(ODEs_DOP853_coe_d614);
    arb_init(ODEs_DOP853_coe_d615);
    arb_init(ODEs_DOP853_coe_d616);
    arb_init(ODEs_DOP853_coe_d701);
    arb_init(ODEs_DOP853_coe_d706);
    arb_init(ODEs_DOP853_coe_d707);
    arb_init(ODEs_DOP853_coe_d708);
    arb_init(ODEs_DOP853_coe_d709);
    arb_init(ODEs_DOP853_coe_d710);
    arb_init(ODEs_DOP853_coe_d711);
    arb_init(ODEs_DOP853_coe_d712);
    arb_init(ODEs_DOP853_coe_d713);
    arb_init(ODEs_DOP853_coe_d714);
    arb_init(ODEs_DOP853_coe_d715);
    arb_init(ODEs_DOP853_coe_d716);
    
    //c7 = convert(T2, 1 // 4)
    arb_one(ODEs_DOP853_coe_c7);
    arb_div_ui(ODEs_DOP853_coe_c7,ODEs_DOP853_coe_c7,4,prec);
    
    //c8 = convert(T2, 4 // 13)
    arb_set_ui(ODEs_DOP853_coe_c8,4);
    arb_div_ui(ODEs_DOP853_coe_c8,ODEs_DOP853_coe_c8,13,prec);
    
    //c9 = convert(T2, 127 // 195)
    arb_set_ui(ODEs_DOP853_coe_c9,127);
    arb_div_ui(ODEs_DOP853_coe_c9,ODEs_DOP853_coe_c9,195,prec);
    
    //c10 = convert(T2, 3 // 5)
    arb_set_ui(ODEs_DOP853_coe_c10,3);
    arb_div_ui(ODEs_DOP853_coe_c10,ODEs_DOP853_coe_c10,5,prec);
    
    //c11 = convert(T2, 6 // 7)
    arb_set_ui(ODEs_DOP853_coe_c11,6);
    arb_div_ui(ODEs_DOP853_coe_c11,ODEs_DOP853_coe_c11,7,prec);
    
    //c6 = convert(T2, 4 // 3 * c7)
    arb_set_ui(ODEs_DOP853_coe_c6,4);
    arb_div_ui(ODEs_DOP853_coe_c6,ODEs_DOP853_coe_c6,3,prec);
    arb_mul(ODEs_DOP853_coe_c6,ODEs_DOP853_coe_c6,ODEs_DOP853_coe_c7,prec);
    
    //c5 = convert(T2, (6 + sqrt(6)) / 10 * c6)
    arb_set_ui(ODEs_DOP853_coe_c5,6);
    arb_sqrt(ODEs_DOP853_coe_c5,ODEs_DOP853_coe_c5,prec);
    arb_add_ui(ODEs_DOP853_coe_c5,ODEs_DOP853_coe_c5,6,prec);
    arb_div_ui(ODEs_DOP853_coe_c5,ODEs_DOP853_coe_c5,10,prec);
    arb_mul(ODEs_DOP853_coe_c5,ODEs_DOP853_coe_c5,ODEs_DOP853_coe_c6,prec);
    
    //c4 = convert(T2, (6 - sqrt(6)) / 10 * c6)
    arb_set_ui(ODEs_DOP853_coe_c4,6);
    arb_sqrt(ODEs_DOP853_coe_c4,ODEs_DOP853_coe_c4,prec);
    arb_neg(ODEs_DOP853_coe_c4,ODEs_DOP853_coe_c4);
    arb_add_ui(ODEs_DOP853_coe_c4,ODEs_DOP853_coe_c4,6,prec);
    arb_div_ui(ODEs_DOP853_coe_c4,ODEs_DOP853_coe_c4,10,prec);
    arb_mul(ODEs_DOP853_coe_c4,ODEs_DOP853_coe_c4,ODEs_DOP853_coe_c6,prec);
    
    //c3 = convert(T2, 2 // 3 * c4)
    arb_set_ui(ODEs_DOP853_coe_c3,2);
    arb_div_ui(ODEs_DOP853_coe_c3,ODEs_DOP853_coe_c3,3,prec);
    arb_mul(ODEs_DOP853_coe_c3,ODEs_DOP853_coe_c3,ODEs_DOP853_coe_c4,prec);
    
    //c2 = convert(T2, 2 // 3 * c3)
    arb_set_ui(ODEs_DOP853_coe_c2,2);
    arb_div_ui(ODEs_DOP853_coe_c2,ODEs_DOP853_coe_c2,3,prec);
    arb_mul(ODEs_DOP853_coe_c2,ODEs_DOP853_coe_c2,ODEs_DOP853_coe_c3,prec);
    
    arb_set_str(ODEs_DOP853_coe_b1, "5.42937341165687622380535766363e-2", prec);
    arb_set_str(ODEs_DOP853_coe_b6, "4.45031289275240888144113950566", prec);
    arb_set_str(ODEs_DOP853_coe_b7, "1.89151789931450038304281599044", prec);
    arb_set_str(ODEs_DOP853_coe_b8, "-5.8012039600105847814672114227", prec);
    arb_set_str(ODEs_DOP853_coe_b9, "3.1116436695781989440891606237e-1", prec);
    arb_set_str(ODEs_DOP853_coe_b10, "-1.52160949662516078556178806805e-1", prec);
    arb_set_str(ODEs_DOP853_coe_b11, "2.01365400804030348374776537501e-1", prec);
    arb_set_str(ODEs_DOP853_coe_b12, "4.47106157277725905176885569043e-2", prec);
    
    arb_set_str(ODEs_DOP853_coe_bhh1, "0.244094488188976377952755905512", prec);
    arb_set_str(ODEs_DOP853_coe_bhh2, "0.733846688281611857341361741547", prec);
    arb_set_str(ODEs_DOP853_coe_bhh3, "0.220588235294117647058823529412e-01", prec);
    
    arb_set_str(ODEs_DOP853_coe_btilde1, "-1.898007540724076157147023288757e-1", prec);
    arb_set_str(ODEs_DOP853_coe_btilde6, "4.45031289275240888144113950566", prec);
    arb_set_str(ODEs_DOP853_coe_btilde7, "1.89151789931450038304281599044", prec);
    arb_set_str(ODEs_DOP853_coe_btilde8, "-5.8012039600105847814672114227", prec);
    arb_set_str(ODEs_DOP853_coe_btilde9, "-4.22682321323791962932445679177e-1", prec);
    arb_set_str(ODEs_DOP853_coe_btilde10, "-1.52160949662516078556178806805e-1", prec);
    arb_set_str(ODEs_DOP853_coe_btilde11, "2.01365400804030348374776537501e-1", prec);
    arb_set_str(ODEs_DOP853_coe_btilde12, "2.26517921983608258118062039631e-2", prec);
    
    arb_set_str(ODEs_DOP853_coe_er1, "0.1312004499419488073250102996e-01", prec);
    arb_set_str(ODEs_DOP853_coe_er6, "-0.1225156446376204440720569753e+01", prec);
    arb_set_str(ODEs_DOP853_coe_er7, "-0.4957589496572501915214079952", prec);
    arb_set_str(ODEs_DOP853_coe_er8, "0.1664377182454986536961530415e+01", prec);
    arb_set_str(ODEs_DOP853_coe_er9, "-0.3503288487499736816886487290", prec);
    arb_set_str(ODEs_DOP853_coe_er10, "0.3341791187130174790297318841", prec);
    arb_set_str(ODEs_DOP853_coe_er11, "0.8192320648511571246570742613e-01", prec);
    arb_set_str(ODEs_DOP853_coe_er12, "-0.2235530786388629525884427845e-01", prec);
    
    arb_set_str(ODEs_DOP853_coe_a0201, "5.26001519587677318785587544488e-2", prec);
    arb_set_str(ODEs_DOP853_coe_a0301, "1.97250569845378994544595329183e-2", prec);
    arb_set_str(ODEs_DOP853_coe_a0302, "5.91751709536136983633785987549e-2", prec);
    arb_set_str(ODEs_DOP853_coe_a0401, "2.95875854768068491816892993775e-2", prec);
    arb_set_str(ODEs_DOP853_coe_a0403, "8.87627564304205475450678981324e-2", prec);
    arb_set_str(ODEs_DOP853_coe_a0501, "2.41365134159266685502369798665e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a0503, "-8.84549479328286085344864962717e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a0504, "9.24834003261792003115737966543e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a0601, "3.7037037037037037037037037037e-2", prec);
    arb_set_str(ODEs_DOP853_coe_a0604, "1.70828608729473871279604482173e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a0605, "1.25467687566822425016691814123e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a0701, "3.7109375e-2", prec);
    arb_set_str(ODEs_DOP853_coe_a0704, "1.70252211019544039314978060272e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a0705, "6.02165389804559606850219397283e-2", prec);
    arb_set_str(ODEs_DOP853_coe_a0706, "-1.7578125e-2", prec);
    
    arb_set_str(ODEs_DOP853_coe_a0801, "3.70920001185047927108779319836e-2", prec);
    arb_set_str(ODEs_DOP853_coe_a0804, "1.70383925712239993810214054705e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a0805, "1.07262030446373284651809199168e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a0806, "-1.53194377486244017527936158236e-2", prec);
    arb_set_str(ODEs_DOP853_coe_a0807, "8.27378916381402288758473766002e-3", prec);
    arb_set_str(ODEs_DOP853_coe_a0901, "6.24110958716075717114429577812e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a0904, "-3.36089262944694129406857109825", prec);
    arb_set_str(ODEs_DOP853_coe_a0905, "-8.68219346841726006818189891453e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a0906, "2.75920996994467083049415600797e1", prec);
    arb_set_str(ODEs_DOP853_coe_a0907, "2.01540675504778934086186788979e1", prec);
    arb_set_str(ODEs_DOP853_coe_a0908, "-4.34898841810699588477366255144e1", prec);
    arb_set_str(ODEs_DOP853_coe_a1001, "4.77662536438264365890433908527e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a1004, "-2.48811461997166764192642586468e0", prec);
    arb_set_str(ODEs_DOP853_coe_a1005, "-5.90290826836842996371446475743e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a1006, "2.12300514481811942347288949897e1", prec);
    arb_set_str(ODEs_DOP853_coe_a1007, "1.52792336328824235832596922938e1", prec);
    arb_set_str(ODEs_DOP853_coe_a1008, "-3.32882109689848629194453265587e1", prec);
    arb_set_str(ODEs_DOP853_coe_a1009, "-2.03312017085086261358222928593e-2", prec);
    
    
    arb_set_str(ODEs_DOP853_coe_a1101, "-9.3714243008598732571704021658e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a1104, "5.18637242884406370830023853209e0", prec);
    arb_set_str(ODEs_DOP853_coe_a1105, "1.09143734899672957818500254654e0", prec);
    arb_set_str(ODEs_DOP853_coe_a1106, "-8.14978701074692612513997267357e0", prec);
    arb_set_str(ODEs_DOP853_coe_a1107, "-1.85200656599969598641566180701e1", prec);
    arb_set_str(ODEs_DOP853_coe_a1108, "2.27394870993505042818970056734e1", prec);
    arb_set_str(ODEs_DOP853_coe_a1109, "2.49360555267965238987089396762e0", prec);
    arb_set_str(ODEs_DOP853_coe_a1110, "-3.0467644718982195003823669022e0", prec);
    arb_set_str(ODEs_DOP853_coe_a1201, "2.27331014751653820792359768449e0", prec);
    arb_set_str(ODEs_DOP853_coe_a1204, "-1.05344954667372501984066689879e1", prec);
    arb_set_str(ODEs_DOP853_coe_a1205, "-2.00087205822486249909675718444e0", prec);
    arb_set_str(ODEs_DOP853_coe_a1206, "-1.79589318631187989172765950534e1", prec);
    arb_set_str(ODEs_DOP853_coe_a1207, "2.79488845294199600508499808837e1", prec);
    arb_set_str(ODEs_DOP853_coe_a1208, "-2.85899827713502369474065508674e0", prec);
    arb_set_str(ODEs_DOP853_coe_a1209, "-8.87285693353062954433549289258e0", prec);
    arb_set_str(ODEs_DOP853_coe_a1210, "1.23605671757943030647266201528e1", prec);
    arb_set_str(ODEs_DOP853_coe_a1211, "6.43392746015763530355970484046e-1", prec);
    
    //c14 = convert(T2, 1 // 10)
    arb_set_ui(ODEs_DOP853_coe_c14,1);
    arb_div_ui(ODEs_DOP853_coe_c14,ODEs_DOP853_coe_c14,10,prec);
    
    //c15 = convert(T2, 2 // 10)
    arb_set_ui(ODEs_DOP853_coe_c15,2);
    arb_div_ui(ODEs_DOP853_coe_c15,ODEs_DOP853_coe_c15,10,prec);
    
    
    //c16 = convert(T2, 7 // 9)
    arb_set_ui(ODEs_DOP853_coe_c16,7);
    arb_div_ui(ODEs_DOP853_coe_c16,ODEs_DOP853_coe_c16,9,prec);
    
    arb_set_str(ODEs_DOP853_coe_a1401, "5.61675022830479523392909219681e-2", prec);
    arb_set_str(ODEs_DOP853_coe_a1407, "2.53500210216624811088794765333e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a1408, "-2.46239037470802489917441475441e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a1409, "-1.24191423263816360469010140626e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a1410, "1.5329179827876569731206322685e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a1411, "8.20105229563468988491666602057e-3", prec);
    arb_set_str(ODEs_DOP853_coe_a1412, "7.56789766054569976138603589584e-3", prec);
    arb_set_str(ODEs_DOP853_coe_a1413, "-8.298e-3", prec);
    
    arb_set_str(ODEs_DOP853_coe_a1501, "3.18346481635021405060768473261e-2", prec);
    arb_set_str(ODEs_DOP853_coe_a1506, "2.83009096723667755288322961402e-2", prec);
    arb_set_str(ODEs_DOP853_coe_a1507, "5.35419883074385676223797384372e-2", prec);
    arb_set_str(ODEs_DOP853_coe_a1508, "-5.49237485713909884646569340306e-2", prec);
    arb_set_str(ODEs_DOP853_coe_a1511, "-1.08347328697249322858509316994e-4", prec);
    arb_set_str(ODEs_DOP853_coe_a1512, "3.82571090835658412954920192323e-4", prec);
    arb_set_str(ODEs_DOP853_coe_a1513, "-3.40465008687404560802977114492e-4", prec);
    arb_set_str(ODEs_DOP853_coe_a1514, "1.41312443674632500278074618366e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a1601, "-4.28896301583791923408573538692e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a1606, "-4.69762141536116384314449447206e0", prec);
    arb_set_str(ODEs_DOP853_coe_a1607, "7.68342119606259904184240953878e0", prec);
    arb_set_str(ODEs_DOP853_coe_a1608, "4.06898981839711007970213554331e0", prec);
    arb_set_str(ODEs_DOP853_coe_a1609, "3.56727187455281109270669543021e-1", prec);
    arb_set_str(ODEs_DOP853_coe_a1613, "-1.39902416515901462129418009734e-3", prec);
    arb_set_str(ODEs_DOP853_coe_a1614, "2.9475147891527723389556272149e0", prec);
    arb_set_str(ODEs_DOP853_coe_a1615, "-9.15095847217987001081870187138e0", prec);
    
    arb_set_str(ODEs_DOP853_coe_d401, "-0.84289382761090128651353491142e+01", prec);
    arb_set_str(ODEs_DOP853_coe_d406, "0.56671495351937776962531783590e+00", prec);
    arb_set_str(ODEs_DOP853_coe_d407, "-0.30689499459498916912797304727e+01", prec);
    arb_set_str(ODEs_DOP853_coe_d408, "0.23846676565120698287728149680e+01", prec);
    arb_set_str(ODEs_DOP853_coe_d409, "0.21170345824450282767155149946e+01", prec);
    arb_set_str(ODEs_DOP853_coe_d410, "-0.87139158377797299206789907490e+00", prec);
    arb_set_str(ODEs_DOP853_coe_d411, "0.22404374302607882758541771650e+01", prec);
    arb_set_str(ODEs_DOP853_coe_d412, "0.63157877876946881815570249290e+00", prec);
    arb_set_str(ODEs_DOP853_coe_d413, "-0.88990336451333310820698117400e-01", prec);
    arb_set_str(ODEs_DOP853_coe_d414, "0.18148505520854727256656404962e+02", prec);
    arb_set_str(ODEs_DOP853_coe_d415, "-0.91946323924783554000451984436e+01", prec);
    arb_set_str(ODEs_DOP853_coe_d416, "-0.44360363875948939664310572000e+01", prec);
    
    
    arb_set_str(ODEs_DOP853_coe_d501, "0.10427508642579134603413151009e+02", prec);
    arb_set_str(ODEs_DOP853_coe_d506, "0.24228349177525818288430175319e+03", prec);
    arb_set_str(ODEs_DOP853_coe_d507, "0.16520045171727028198505394887e+03", prec);
    arb_set_str(ODEs_DOP853_coe_d508, "-0.37454675472269020279518312152e+03", prec);
    arb_set_str(ODEs_DOP853_coe_d509, "-0.22113666853125306036270938578e+02", prec);
    arb_set_str(ODEs_DOP853_coe_d510, "0.77334326684722638389603898808e+01", prec);
    arb_set_str(ODEs_DOP853_coe_d511, "-0.30674084731089398182061213626e+02", prec);
    arb_set_str(ODEs_DOP853_coe_d512, "-0.93321305264302278729567221706e+01", prec);
    arb_set_str(ODEs_DOP853_coe_d513, "0.15697238121770843886131091075e+02", prec);
    arb_set_str(ODEs_DOP853_coe_d514, "-0.31139403219565177677282850411e+02", prec);
    arb_set_str(ODEs_DOP853_coe_d515, "-0.93529243588444783865713862664e+01", prec);
    arb_set_str(ODEs_DOP853_coe_d516, "0.35816841486394083752465898540e+02", prec);
    
    arb_set_str(ODEs_DOP853_coe_d601, "0.19985053242002433820987653617e+02", prec);
    arb_set_str(ODEs_DOP853_coe_d606, "-0.38703730874935176555105901742e+03", prec);
    arb_set_str(ODEs_DOP853_coe_d607, "-0.18917813819516756882830838328e+03", prec);
    arb_set_str(ODEs_DOP853_coe_d608, "0.52780815920542364900561016686e+03", prec);
    arb_set_str(ODEs_DOP853_coe_d609, "-0.11573902539959630126141871134e+02", prec);
    arb_set_str(ODEs_DOP853_coe_d610, "0.68812326946963000169666922661e+01", prec);
    arb_set_str(ODEs_DOP853_coe_d611, "-0.10006050966910838403183860980e+01", prec);
    arb_set_str(ODEs_DOP853_coe_d612, "0.77771377980534432092869265740e+00", prec);
    arb_set_str(ODEs_DOP853_coe_d613, "-0.27782057523535084065932004339e+01", prec);
    arb_set_str(ODEs_DOP853_coe_d614, "-0.60196695231264120758267380846e+02", prec);
    arb_set_str(ODEs_DOP853_coe_d615, "0.84320405506677161018159903784e+02", prec);
    arb_set_str(ODEs_DOP853_coe_d616, "0.11992291136182789328035130030e+02", prec);
    
    arb_set_str(ODEs_DOP853_coe_d701, "-0.25693933462703749003312586129e+02", prec);
    arb_set_str(ODEs_DOP853_coe_d706, "-0.15418974869023643374053993627e+03", prec);
    arb_set_str(ODEs_DOP853_coe_d707, "-0.23152937917604549567536039109e+03", prec);
    arb_set_str(ODEs_DOP853_coe_d708, "0.35763911791061412378285349910e+03", prec);
    arb_set_str(ODEs_DOP853_coe_d709, "0.93405324183624310003907691704e+02", prec);
    arb_set_str(ODEs_DOP853_coe_d710, "-0.37458323136451633156875139351e+02", prec);
    arb_set_str(ODEs_DOP853_coe_d711, "0.10409964950896230045147246184e+03", prec);
    arb_set_str(ODEs_DOP853_coe_d712, "0.29840293426660503123344363579e+02", prec);
    arb_set_str(ODEs_DOP853_coe_d713, "-0.43533456590011143754432175058e+02", prec);
    arb_set_str(ODEs_DOP853_coe_d714, "0.96324553959188282948394950600e+02", prec);
    arb_set_str(ODEs_DOP853_coe_d715, "-0.39177261675615439165231486172e+02", prec);
    arb_set_str(ODEs_DOP853_coe_d716, "-0.14972683625798562581422125276e+03", prec);
    
    ODEs_get_DOP853_cal_coe_help=_arb_vec_init(1); //用于判断初始化是否完成
    arb_one(ODEs_get_DOP853_cal_coe_help);
}
