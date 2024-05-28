#ifndef __PBHS_OTHER_OUTPUT__   /* Include guard */
#define __PBHS_OTHER_OUTPUT__

#include <arb.h> //高精度实数运算
#include "../new_type.h"
#include "../phy_constant.h"

int get_save_path(char* res); //获得存储文件路径
void print_progress(size_t count, size_t max); //进度条显示


#endif // __PBHS_OTHER_OUTPUT__ 
