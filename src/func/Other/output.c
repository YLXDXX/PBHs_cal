#include "output.h" 
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <string.h> 
#include <stdlib.h>


//获取存储文件路径
int get_save_path(char* res)
{
    //获取程序运行路径，对其进行处理后，用来确定存储一些文件的路径
    
    
    char run_cwd[PATH_MAX+1]; // 其中 PATH_MAX 路径名最大字节数 NAME_MAX 文件名最大字节数
    
    if ( readlink("/proc/self/exe", run_cwd, sizeof(run_cwd) ) != 0) //通过 /proc/self/exe 获取程序的名称和绝对路径
                                                                    //注意，通过 getcwd() 获取的结果跟执行情况有关
    {
        if(Stdout_verbose==true)
        {
            printf("\nRunning command dir: %s\n", run_cwd);
        }
        
    } else
    {
        perror("readlink() error");
        exit(1);
    }
    
    // 路径示例 /home/shui/Downloads/PBHs/build/linux/x86_64/release/PBHs
    char delim[] = "/";  //分割符
    char* save_split[PATH_MAX+1];//此数组用于存分割后的内容
    char* temp_buffer = strdup(run_cwd); //给 strsep 函数用
    char* token; 
    int count=0,i=0;
    
    while((token = strsep(&temp_buffer, delim)) && count < PATH_MAX+1)
    {
        save_split[count++] = token;
    }
    
    free(temp_buffer);
    
    
    //这里路径的取法是相对于xmake来的
    //char res[PATH_MAX+1]=""; //获取最终路径
    for (i = 0; i < count-5; i++)  //当前目录需要向上走 5 级，得 /home/shui/Downloads/PBHs
    {
        if( ! strcmp(save_split[i], "") ) //字符串比较，相同返回 0 
        {
            strcat(res, delim); //字符串添加
            
        }else
        {
            strcat(res, save_split[i]); //字符串添加
            strcat(res, delim);
        }
    }
    
    strcpy(Out_date_file, Path_save); //字符串复制
    strcpy(Out_fitted_file, Path_save);
    strcpy(Out_fitted_x, Path_save);
    strcpy(Out_fitted_y, Path_save);
    strcpy(Out_picture_file, Path_save);
    
    
    if(Stdout_verbose==true)
    {
        printf("Current saving dir: %s\n\n", res);
    }
    
    return 0; 
}



//输出进度条显示
//见 https://gist.github.com/amullins83/24b5ef48657c08c4005a8fab837b7499
void print_progress(size_t count, size_t max) {
    const int bar_width = 50;
    
    float progress = (float) count / max;
    int bar_length = progress * bar_width;
    
    printf("\rProgress: [");
    for (int i = 0; i < bar_length; ++i) {
        printf("#");
    }
    for (int i = bar_length; i < bar_width; ++i) {
        printf(" ");
    }
    printf("] %.2f%%", progress * 100);
    
    fflush(stdout);
}


//将一组矢量 x,y 中每点的值输出到文件
void Vector_point_output_to_file(const arb_ptr x, const arb_ptr y, const slong num, char tag)
{
    char* mod=&tag;
    int digit=18; //输出数字的有效位数
    
    FILE * fp;
    fp = fopen(Out_picture_file, mod); //打开文件，a追加，w重新写入
    
    if( fp == NULL ) { //对文件打开操作进行判断
        printf("\n\nOpen Error: %s\t\n",Out_picture_file);perror("file");printf("\n");
        exit(-1);
    }
    
    for(slong i=0; i<num; i++)
    {
        arb_fprintn(fp,x+i,digit,ARB_STR_NO_RADIUS);
        fprintf(fp, "\t");
        arb_fprintn(fp,y+i,digit,ARB_STR_NO_RADIUS);
        fprintf(fp, "\n");
    }
    
    fclose(fp); //关闭文件
    printf("\n输出完成，共 %li 个点 \n",num);
    
}

//将一组矢量 x,y 中每点的值输出到文件，输出格式保持arb内部的表达式
void Vector_point_write_to_file_arb(const arb_ptr x, const arb_ptr y, const slong num, slong prec)
{
    arb_t ln_x,ln_y;
    arb_init(ln_x);
    arb_init(ln_y);
    
    FILE * fp_x;
    FILE * fp_y;
    
    fp_x = fopen(Out_fitted_x, "w"); //打开文件，a追加，w重新写入
    fp_y = fopen(Out_fitted_y, "w");
    
    if( fp_x == NULL || fp_y == NULL ) { //对文件打开操作进行判断
        printf("\n\nOpen Error: %s\t\n",Out_fitted_y);perror("file");printf("\n");
        exit(-1);
    }
    
    for(slong i=0; i<num; i++)
    {
        arb_log(ln_x,x+i,prec); //对 x 取对数 ln(x)
        arb_dump_file(fp_x,ln_x);
        fprintf(fp_x, "\n");
        
        //arb_log(ln_y,y+i,prec);
        arb_dump_file(fp_y,y+i);
        fprintf(fp_y, "\n");
    }
    
    arb_clear(ln_x);
    arb_clear(ln_y);
    
    fclose(fp_x); //关闭文件
    fclose(fp_y);
    
    printf("\n写入完成，共 %li 个点 \n",num);
    
}


void Vector_point_read_to_file_arb(const arb_ptr x, const arb_ptr y, const slong num, slong prec)
{
    FILE * fp_x;
    FILE * fp_y;
    
    fp_x = fopen(Out_fitted_x, "r"); //打开文件，a追加，w重新写入
    fp_y = fopen(Out_fitted_y, "r");
    
    if( fp_x == NULL || fp_y == NULL ) { //对文件打开操作进行判断
        printf("\n\nOpen Error: %s\t\n",Out_fitted_y);perror("file");printf("\n");
        exit(-1);
    }
    
    for(slong i=0; i<num; i++)
    {
        arb_load_file(x+i,fp_x);
        
        arb_load_file(y+i,fp_y);
    }
    
    fclose(fp_x); //关闭文件
    fclose(fp_y);
    
    printf("\n读取完成，共 %li 个点 \n",num);
    
}
