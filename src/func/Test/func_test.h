#ifndef __PBHS_TEST_FUNC_TEST__   /* Include guard */
#define __PBHS_TEST_FUNC_TEST__

#include <arb.h> //高精度实数运算

#include "../new_type.h"
#include "../phy_constant.h"
#include "../math_method.h"



int Func_test(arb_t res, const arb_t x, void* params, const slong order, slong prec); //测试函数，用于各种测试

int Func_test_2(arb_t res, const arb_t x, const arb_t y, void* params, const slong order, slong prec); //测试函数，用于各种测试

int Func_test_3(arb_t res, const arb_t p1, const arb_t p2, void* params, const slong order, slong prec); //测试函数，用于各种测试

int Find_root_test(arb_t res, const arb_t x, void* params, const slong order, slong prec);

int Func_test_4(arb_t res, const arb_t x, void* params, const slong order, slong prec); //测试函数，用于各种测试

int Func_test_5(arb_t res, const arb_t x, void* params, const slong order, slong prec); //测试函数，用于各种测试

int Func_test_6(arb_t res, const arb_t x, void* params, const slong order, slong prec); //测试函数，用于各种测试

int Func_test_7(arb_t res, const arb_t x, void* params, const slong order, slong prec); //测试函数，用于各种测试

int Func_test_8(arb_t res, const arb_t x, void* params, const slong order, slong prec);

int Func_test_9(arb_t res, const arb_t x, void* params, const slong order, slong prec);

int Func_test_10(arb_t res, const arb_t x, void* params, const slong order, slong prec);

int Func_test_quad_func_01(arb_t res, const arb_t x, const arb_t y, void* params, const slong order, slong prec);
int Func_test_quad_func_01_y_a(arb_t res, const arb_t x, void* params, const slong order, slong prec);
int Func_test_quad_func_01_y_b(arb_t res, const arb_t x, void* params, const slong order, slong prec);

int Func_test_quad_func_02(arb_t res, const arb_t x, const arb_t y, void* params, const slong order, slong prec);
int Func_test_quad_func_02_y_a(arb_t res, const arb_t x, void* params, const slong order, slong prec);
int Func_test_quad_func_02_y_b(arb_t res, const arb_t x, void* params, const slong order, slong prec);


int Func_test_quad_func_03(arb_t res, const arb_t x, const arb_t y, void* params, const slong order, slong prec);
int Func_test_quad_func_03_y_a(arb_t res, const arb_t x, void* params, const slong order, slong prec);
int Func_test_quad_func_03_y_b(arb_t res, const arb_t x, void* params, const slong order, slong prec);



#endif // __PBHS_TEST_FUNC_TEST__  
