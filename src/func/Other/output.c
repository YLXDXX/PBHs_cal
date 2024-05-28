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


