#ifndef __PBHS_OTHER_OUTPUT__   /* Include guard */
#define __PBHS_OTHER_OUTPUT__

#include <arb.h> //高精度实数运算
#include "../new_type.h"
#include "../phy_constant.h"

int get_save_path(char* res); //获得存储文件路径
void print_progress(size_t count, size_t max); //进度条显示

//将一组矢量 x,y 中每点的值输出到文件
void Vector_point_output_to_file(const arb_ptr x, const arb_ptr y, const slong num, char tag);

void Vector_point_write_to_file_arb(const arb_ptr x, const arb_ptr y, const slong num,  const int i, slong prec);
void Vector_point_read_to_file_arb(const arb_ptr x, const arb_ptr y, const slong num, const int i, slong prec);

#endif // __PBHS_OTHER_OUTPUT__ 
