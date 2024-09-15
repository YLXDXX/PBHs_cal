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



