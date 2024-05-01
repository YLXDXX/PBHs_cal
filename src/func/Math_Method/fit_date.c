#include "fit_date.h"
#include <arb_mat.h>
#include <profiler.h> //flint自带的多线程
#include <stdlib.h>
#include <string.h>


//函数对y=f(x) 将点 (x,y) 输出到文件，储存备用
int Func_output_point(int (*func)(arb_t f_res, const arb_t f_x, const slong f_order, slong prec),
                      const slong order,
                      const arb_t a, const arb_t b, const slong N, char* file, slong prec)
{
    arb_t x,y,gap;
    arb_init(x);
    arb_init(y);
    arb_init(gap);
    
    arb_sub(gap,b,a,prec); //获得区间间隔
    arb_div_ui(gap,gap,N,prec); 
    
    //计算 (x,y), 然后存储
    slong i;
    
    FILE* fp;
    fp = fopen(file, "w");
    
    if (!fp) {
        printf("Math_Method -> fit_date -> ERROR: Unable to open file: %s", file);
        exit(1);
        return EXIT_FAILURE;
    }
    
    for(i=0; i <= N; i++ ) //倒序输出，读入还会倒一下，就正了
    {
        arb_mul_si(x,gap,i,prec); //更新 x
        arb_sub(x,b,x,prec);
        
        func(y,x,order,prec); //得到函数值
        
        //存储，一个 (x,y) 点一行，用 tab 分隔
        arb_dump_file(fp,x); 
        fprintf(fp, "\t");
        arb_dump_file(fp,y);
        fprintf(fp, "\n");
        
        printf("%li/%li\n",i,N);
    }
    
    fclose(fp);
    
    printf("计算数据写入成功：%s\n\n", file);
    
    arb_clear(x);
    arb_clear(y);
    arb_clear(gap);
    
    return EXIT_SUCCESS;
}


//文件读取，数据拟合
int Func_fit_point(char* fitted_file, char* output_file,const unsigned int num_threads, slong prec)
{
    //https://linux.cn/article-11950-1.html
    //https://opensource.com/article/20/2/c-data-science
    
    //定义一个由 FUNC_INPUT_POINT 构成的单链列表
    SLIST_HEAD(func_input_list, FUNC_INPUT_POINT) head = SLIST_HEAD_INITIALIZER(head);
    SLIST_INIT(&head);
    
    printf("#### 开始读取数据 ####\n");
    
    //读取还应的数据文件
    FILE* input_file = fopen(output_file, "r"); //打开文件
    
    if (!input_file) {
        printf("ERROR: Unable to open file: %s", output_file);
        
        return EXIT_FAILURE;
    }
    
    ulong row = 0; //有效行数，计算成功读取的行数，里面可能有空行等情况
    ulong row_s = 0; //总行数，用以把最开始的几行注释忽略
    
    char* str_x="";
    char* str_y="";
    
    //然后逐行读取文件，直到出现错误或文件结束
    while (!ferror(input_file) && !feof(input_file))
    {
        size_t buffer_size = 0; //size_t是无符号的长整型
                                //一般用来表示字节数的多少。常用于如sizeof返回值的类型
                                //size_t是一种跟具体的平台有关联的类型，会具体调整其能表示的范围，因此其可移植性会更好
                                //size_t能保证可以存储任何类型理论上可能的对象的最大值，包括数组类型。
        char *buffer = NULL;
        
        //getline() 函数可读取文件中的整行，并负责分配必要的内存
        getline(&buffer, &buffer_size, input_file);
        
        
        //读取整行后，使用 strtok() 函数将每一行分成字元「token」。遍历字元，选择所需的列
        //这里理应 strlen(buffer) > 0 但最后一行会遇到一个奇怪的字符，导至多出一行
        //现暴力的将其设为 > 6 个字符长度 
        // 空行 buffer[0] !='\n' ，用前面6个字符长度已不需要
        if (row_s >= Func_output_skip_header && strlen(buffer) > 6)
        {
            unsigned int column = 0; //对于列数来说，足够
            
            char* token = strtok(buffer, Func_output_delimiter);
            
            while (token != NULL) //读取一行中的各例，读取的内容是字符串的形式
                //可通过 sscanf(token, "%lf", &value); 的形式
                //把字符串的内容给 double value 变量
            {
                //这里输入的是字符串
                //由 int arb_load_str(arb_t x, const char *str) 直接读入即可
                if (column == Func_output_column_x) {
                    str_x = token;
                } else if (column == Func_output_column_y) {
                    str_y = token;
                }
                
                column += 1;
                token = strtok(NULL, Func_output_delimiter);
            }
            
            //当选择了 x 和 y 值时，将新数据点插入链表中：
            //malloc() 函数为新数据点动态分配（保留）一些持久性内存
            //struct FUNC_INPUT_POINT *datum = malloc(sizeof(struct FUNC_INPUT_POINT));
            struct FUNC_INPUT_POINT *datum = (struct FUNC_INPUT_POINT *)calloc(1,sizeof(struct FUNC_INPUT_POINT));
            
            int jude_1,jude_2;
            
            //arb_t 使用时，需初始化，即使是结构体，也得初始化
            arb_init(datum->x); //arb_t 使用前初始化，避免可能存在的内存分配问题
                                //如没初始化时，当其遇是整数的时候，不会按 arb_t 结构体 48 分配空间
                                //这样，就导制存存空间分配的问题，发生段错误 Segmentation fault (core dumped)
            arb_init(datum->y);
            
            jude_1=arb_load_str(datum->x,str_x);
            jude_2=arb_load_str(datum->y,str_y);
            
            if( jude_1 != 0 || jude_2 != 0)
            {
                printf("Math_Method -> fit_date -> 加入数据失败，行数：%li，文件名：%s\n",row_s,Func_output_file);
                exit(1);
            }
            
            SLIST_INSERT_HEAD(&head, datum, func_input_point_entries);
            
            row+= 1; //有效行数
        }
        
        //释放内存，改变行数，读取下一行
        free(buffer);
        row_s += 1; //总行数
    }
    
    fclose(input_file); //读取文件完成，关闭文件
    
    if(row==0)
    {
        printf("读取错误，条目为 0 ：%s\n",output_file);
        return EXIT_FAILURE;
    }
    
    
    //由于你将不知道要创建的数组的大小，因此必须手动分配它们的内存
    //注意，这里分配的内存，即数组中变量的个数，需与数据的个数严格匹配
    
    const size_t entries_number = row; //计算所有数据行数，得到娈量个数
    
    //将相应的数据用存到数组中，供下一步使用
    //手动分配内存的方法有两种 calloc 和 malloc ，malloc在数据点很少「如10」时有问题
    arb_t* array_x =(arb_t *) calloc(entries_number,sizeof(arb_t)); //根据变量的数目分配内存
    arb_t* array_y =(arb_t *) calloc(entries_number,sizeof(arb_t));
    
    //arb_t* array_x =(arb_t *) malloc(entries_number * sizeof(arb_t)); //根据变量的数目分配内存
    //arb_t* array_y =(arb_t *) malloc(entries_number *sizeof(arb_t));
    
    
    //arb_t 使用时，需初始化，即使是数组，也得初始化
    for(ulong i=0; i < entries_number; i+=1)
    {
        arb_init(array_x[i]);
        arb_init(array_y[i]);
    }
    
    
    
    if (!array_x || !array_y) {
        printf("ERROR: Unable to allocate data arrays\n");
        return EXIT_FAILURE;
    }
    
    
    struct FUNC_INPUT_POINT *datum;
    datum = SLIST_FIRST(&head);
    
    ulong i = 0; //给数组 x y 赋值用
    
    //遍历链表，存入数组
    SLIST_FOREACH(datum, &head, func_input_point_entries)
    {
        arb_set(array_x[i],datum->x); //存值到数组
        arb_set(array_y[i],datum->y);
        
        i += 1;
    }
    
    //处理完链表后，需清理
    //手动释放已分配的内存，防止内存泄漏
    while (!SLIST_EMPTY(&head))
    {
        struct FUNC_INPUT_POINT *datum = SLIST_FIRST(&head);
        
        SLIST_REMOVE_HEAD(&head, func_input_point_entries);
        
        free(datum);
    }
    
    
    printf("\n%s 读取成功\n有效数据对 (x,y) 共计：%li\n\n", Func_output_file,entries_number);
    
    
    printf("运用三次样条插值，开始求解相关矩阵系数\n\n");
    
    //拟合数据，利用三次三插的办法
    //最后会解线性方程组，一个系数矩阵，两个矢量， AM=B 求解 M 
    
    arb_mat_t M,B; //两个矢量 AM=B
    arb_mat_t A; //系数矩阵
    
    arb_mat_init(M,entries_number-2,1); //由于使用自然边界条件 M需要求的值由 n 个变为 n-2 个
    arb_mat_init(B,entries_number-2,1);
    arb_mat_init(A,entries_number-2,entries_number-2);
    
    //先将所有矩阵系数设为0
    arb_mat_zero(A);
    arb_mat_zero(B);
    arb_mat_zero(M);
    
    
    
    //已知 y_i=f(x_i) i=0,1,2,..,n 给出 n+1 个点
    //先求s_h, s_lambda, s_nu, s_d,这些中间值，再求A、B 两个矩阵
    //这里利用的边界条任是 己知 S''(x_0) = M_0=0, S''(x_n) = M_n=0 「自然边界条件」
    //还需要求M_1，M_2，...，M_{n-1}
    
    arb_t s_h,s_h_sub,s_lambda,s_nu,s_d,d_q_2,d_q_3,temp,temp_2;
    arb_init(s_h);
    arb_init(s_h_sub);
    arb_init(s_lambda);
    arb_init(s_nu);
    arb_init(s_d);
    arb_init(d_q_2);
    arb_init(d_q_3);
    arb_init(temp);
    arb_init(temp_2);
    
    
    for(ulong i=1; i < entries_number-1; i+=1 ) //这时 entries_number-1 和  i=1是由于自然边界条件
    {
        
        
        //求h_i
        arb_sub(s_h,array_x[i+1],array_x[i],prec);
        //h_{i-1}
        arb_sub(s_h_sub,array_x[i],array_x[i-1],prec);
        
        
        //利用h_i 求 lambda_i 
        arb_add(s_lambda,s_h_sub,s_h,prec);
        arb_div(s_lambda,s_h,s_lambda,prec);
        
        //利用lambda_i 求 nu_i
        arb_neg(s_nu,s_lambda);
        arb_add_si(s_nu,s_nu,1,prec);
        
        //求d_i，这时需要计算三阶差商
        //f[x_{i-1}, x_i, x_{i+1}] =f[x_i,x_{i+1}] - f[x_[i-1],x_i] )/(x_{i+1}-x_{i-1})
        //计算二阶差商
        arb_sub(temp,array_y[i+1],array_y[i],prec);
        arb_sub(d_q_2,array_x[i+1],array_x[i],prec);
        arb_div(d_q_2,temp,d_q_2,prec);
        
        arb_sub(temp,array_y[i],array_y[i-1],prec);
        arb_sub(temp_2,array_x[i],array_x[i-1],prec);
        arb_div(temp,temp,temp_2,prec);
        
        //计算三阶差商
        arb_sub(temp,d_q_2,temp,prec);
        arb_sub(temp_2,array_x[i+1],array_x[i-1],prec);
        arb_div(d_q_3,temp,temp_2,prec);
        
        //计算d_i
        arb_div_ui(s_d,d_q_3,6,prec);
        
        
        //辅助变量计算完成，计算相应的矩阵
        //矩阵A
        if (i==1)
        {
            //arb_set(arb_mat_entry(A, i-1, i-2),s_nu);
            arb_set_si(arb_mat_entry(A, i-1, i-1),2);
            arb_set(arb_mat_entry(A, i-1, i),s_lambda);
        }else if (i==entries_number-2)
        {
            arb_set(arb_mat_entry(A, i-1, i-2),s_nu);
            arb_set_si(arb_mat_entry(A, i-1, i-1),2);
            //arb_set(arb_mat_entry(A, i-1, i),s_lambda);
        }else
        {
            arb_set(arb_mat_entry(A, i-1, i-2),s_nu);
            arb_set_si(arb_mat_entry(A, i-1, i-1),2);
            arb_set(arb_mat_entry(A, i-1, i),s_lambda);
            
        }
        
        
        //矩阵B
        arb_set(arb_mat_entry(B, i-1, 0),s_d);
        
    }
    
    arb_clear(s_h);//释放
    arb_clear(s_h_sub);
    arb_clear(s_lambda);
    arb_clear(s_nu);
    arb_clear(s_d);
    arb_clear(d_q_2);
    arb_clear(d_q_3);
    arb_clear(temp);
    arb_clear(temp_2);
    
    printf("相关矩阵系数求解完成\n\n");
    
    printf("开始求解线性方程组\n\n");
    
    flint_set_num_threads(num_threads); //多线程设定
    
    TIMEIT_ONCE_START //打印求解矩阵用时
    
    //此矩阵，适宜用 LU 分解来求解，特别是当点数量比较多时
    arb_mat_solve_lu(M,A,B,prec);
    
    TIMEIT_ONCE_STOP //多线程结束
    
    printf("线性方程组求解完成\n\n");
    
    //arb_printd(arb_mat_entry(B,9,0),5);printf("\n");
    //printf("\n %li X %li\n",arb_mat_nrows(B),arb_mat_ncols(B));
    //arb_mat_printd(D,3);
    
    printf("求解各个区间的三次多项式系数\n\n");
    
    //定义最终存储的结构体数组 arry_fit_save 
    struct FUNC_FITTED_DATE *arry_fit_save =(struct FUNC_FITTED_DATE *) calloc(entries_number,sizeof(struct FUNC_FITTED_DATE)); //根据变量的数目分配内存
    
    
    arb_t s,t,w,h_i;
    
    arb_init(s);
    arb_init(t);
    arb_init(w);
    arb_init(h_i);
    
    for(ulong i=0; i < entries_number; i+=1)
    {
        //arb_t 使用时，需初始化
        arb_init(arry_fit_save[i].a);
        arb_init(arry_fit_save[i].b);
        arb_init(arry_fit_save[i].c);
        arb_init(arry_fit_save[i].d);
        arb_init(arry_fit_save[i].x);
        
        
        //存对应的坐标
        arb_set(arry_fit_save[i].x, array_x[i]);
        
        //对的M而言，由于采用自然边界条件，导致M_0=M_n-1=0，M矩阵中无这两个点
        // M --> M_1, M_2, M_3， ... M_n-3, M_n-2
        slong j=i-1; //用于恢愎 M 的指标
        
        
        
        if (i==0){
            //此时M_0=0;
            
            //求h_i
            arb_sub(h_i,array_x[i+1],array_x[i],prec);
            
            //常数项
            
            //括号前 0
            
            //括号 中
            arb_mul_ui(t,array_y[i],6,prec);
            arb_mul(s,t,array_x[i+1],prec);
            
            //括号 中
            arb_sqr(w,array_x[i],prec); //有个三次方
            arb_mul(w,w,array_x[i],prec);
            arb_mul(w,w,arb_mat_entry(M,j+1,0),prec);
            arb_sub(s,s,w,prec);
            
            //括号后
            arb_sqr(w,h_i,prec);
            arb_mul(w,w,arb_mat_entry(M,j+1,0),prec);
            arb_mul_ui(t,array_y[i+1],6,prec);
            arb_sub(w,w,t,prec);
            arb_mul(w,w,array_x[i],prec);
            arb_add(s,s,w,prec);
            
            arb_mul_ui(w,h_i,6,prec); //最后相除
            arb_div(arry_fit_save[i].d,s,w,prec);
            
            //一次项
            
            //括号前 0
            
            //括号中
            arb_sqr(w,array_x[i],prec);
            arb_mul(w,w,arb_mat_entry(M,j+1,0),prec);
            arb_mul_ui(s,w,3,prec);
            arb_neg(s,s);
            
            //括号中
            arb_mul(w,arb_mat_entry(M,j+1,0),h_i,prec);
            arb_mul(w,w,h_i,prec);
            arb_add(s,s,w,prec);
            
            //括号后
            arb_sub(w,array_y[i],array_y[i+1],prec);
            arb_mul_ui(w,w,6,prec);
            arb_add(s,s,w,prec);
            
            arb_mul_ui(w,h_i,6,prec); //最后相除
            arb_div(arry_fit_save[i].c,s,w,prec);
            arb_neg(arry_fit_save[i].c,arry_fit_save[i].c);
            
            //二次项
            
            arb_mul(w,arb_mat_entry(M,j+1,0),array_x[i],prec);
            arb_neg(s,w),
            
            arb_mul_ui(w,h_i,2,prec); //最后相除
            arb_div(arry_fit_save[i].b,s,w,prec);
            
            //三次项
            arb_set(s,arb_mat_entry(M,j+1,0));
            
            arb_mul_ui(w,h_i,6,prec); //最后相除
            arb_div(arry_fit_save[i].a,s,w,prec);
            
            continue;
        }
        
        if (i==entries_number-2){
            //最后一个插值 [x_n-2, x_n-1] 其n个点，从零开始计数
            // M_n-1=0
            
            //常数项
            
            //括号前
            arb_sqr(s,array_x[i+1],prec); //有个三次方
            arb_mul(s,s,array_x[i+1],prec);
            arb_mul(s,s,arb_mat_entry(M,j,0),prec);
            
            //括号 中
            arb_mul_ui(t,array_y[i],6,prec);
            arb_sqr(w,h_i,prec);
            arb_mul(w,arb_mat_entry(M,j,0),w,prec);
            arb_sub(t,t,w,prec);
            arb_mul(t,t,array_x[i+1],prec);
            arb_add(s,s,t,prec);
            
            //括号后
            arb_mul_si(w,array_y[i+1],-6,prec);
            arb_mul(w,w,array_x[i],prec);
            arb_add(s,s,w,prec);
            
            arb_mul_ui(w,h_i,6,prec); //最后相除
            arb_div(arry_fit_save[i].d,s,w,prec);
            
            //一次项
            
            //括号前
            arb_sqr(s,array_x[i+1],prec);
            arb_mul(s,s,arb_mat_entry(M,j,0),prec);
            arb_mul_ui(s,s,3,prec);
            
            //括号中
            arb_neg(w,arb_mat_entry(M,j,0));
            arb_mul(w,w,h_i,prec);
            arb_mul(w,w,h_i,prec);
            arb_add(s,s,w,prec);
            
            //括号后
            arb_sub(w,array_y[i],array_y[i+1],prec);
            arb_mul_ui(w,w,6,prec);
            arb_add(s,s,w,prec);
            
            arb_mul_ui(w,h_i,6,prec); //最后相除
            arb_div(arry_fit_save[i].c,s,w,prec);
            arb_neg(arry_fit_save[i].c,arry_fit_save[i].c); //最后要取负号
            
            //二次项
            
            arb_mul(s,arb_mat_entry(M,j,0),array_x[i+1],prec);
            
            arb_mul_ui(w,h_i,2,prec); //最后相除
            arb_div(arry_fit_save[i].b,s,w,prec);
            
            //三次项
            
            arb_neg(s, arb_mat_entry(M,j,0));
            
            arb_mul_ui(w,h_i,6,prec); //最后相除
            arb_div(arry_fit_save[i].a,s,w,prec);
            
            continue;
        }
        
        if (i==entries_number-1){
            //最大的点，向后无插值
            continue;
        }
        
        
        
        
        
        //求h_i
        arb_sub(h_i,array_x[i+1],array_x[i],prec);
        
        //求解原函数，一阶导，二阶导
        //S_i: M_i * (x_i1 - x)^3 / (6*h_i) + M_i1 * (x-x_i)^3/(6*h_i) + (f_i - M_i * (h_i)^2/6)*(x_i1-x)/h_i + (f_i1 - M_i1 * (h_i)^2/6)*(x-x_i)/h_i ;
        //maxima 求得 x 的各项系数 coeff(expand(S_i),x,0); coeff(expand(S_i),x,1);
        //
        //另外，此处的 S_i(x) 是小区间 [x_i, x_i+1] 上的三次多项式
        //
        //常数项 (M_i*x_i1^3)/(6*h_i)-(M_i*h_i*x_i1)/6+(f_i*x_i1)/h_i-(M_i1*x_i^3)/(6*h_i)+(M_i1*h_i*x_i)/6-(f_i1*x_i)/h_i
        //     =(M_i*x_i1^3+(6*f_i-M_i*h_i^2)*x_i1-M_i1*x_i^3+(M_i1*h_i^2-6*f_i1)*x_i)/(6*h_i)
        
        
        
        //括号前
        arb_sqr(s,array_x[i+1],prec); //有个三次方
        arb_mul(s,s,array_x[i+1],prec);
        arb_mul(s,s,arb_mat_entry(M,j,0),prec);
        
        //括号 中
        arb_mul_ui(t,array_y[i],6,prec);
        arb_sqr(w,h_i,prec);
        arb_mul(w,arb_mat_entry(M,j,0),w,prec);
        arb_sub(t,t,w,prec);
        arb_mul(t,t,array_x[i+1],prec);
        arb_add(s,s,t,prec);
        
        //括号 中
        arb_sqr(w,array_x[i],prec); //有个三次方
        arb_mul(w,w,array_x[i],prec);
        arb_mul(w,w,arb_mat_entry(M,j+1,0),prec);
        arb_sub(s,s,w,prec);
        
        //括号后
        arb_sqr(w,h_i,prec);
        arb_mul(w,w,arb_mat_entry(M,j+1,0),prec);
        arb_mul_ui(t,array_y[i+1],6,prec);
        arb_sub(w,w,t,prec);
        arb_mul(w,w,array_x[i],prec);
        arb_add(s,s,w,prec);
        
        arb_mul_ui(w,h_i,6,prec); //最后相除
        arb_div(arry_fit_save[i].d,s,w,prec);
        
        
        //一次项 -(M_i*x_i1^2)/(2*h_i)+(M_i1*x_i^2)/(2*h_i)-(M_i1*h_i)/6+(M_i*h_i)/6+f_i1/h_i-f_i/h_i
        //     =-(3*M_i*x_i1^2-3*M_i1*x_i^2+(M_i1-M_i)*h_i^2-6*f_i1+6*f_i)/(6*h_i)
        
        //括号前
        arb_sqr(s,array_x[i+1],prec);
        arb_mul(s,s,arb_mat_entry(M,j,0),prec);
        arb_mul_ui(s,s,3,prec);
        
        //括号中
        arb_sqr(w,array_x[i],prec);
        arb_mul(w,w,arb_mat_entry(M,j+1,0),prec);
        arb_mul_ui(w,w,3,prec);
        arb_sub(s,s,w,prec);
        
        //括号中
        arb_sub(w,arb_mat_entry(M,j+1,0),arb_mat_entry(M,j,0),prec);
        arb_mul(w,w,h_i,prec);
        arb_mul(w,w,h_i,prec);
        arb_add(s,s,w,prec);
        
        //括号后
        arb_sub(w,array_y[i],array_y[i+1],prec);
        arb_mul_ui(w,w,6,prec);
        arb_add(s,s,w,prec);
        
        arb_mul_ui(w,h_i,6,prec); //最后相除
        arb_div(arry_fit_save[i].c,s,w,prec);
        arb_neg(arry_fit_save[i].c,arry_fit_save[i].c); //最后要取负号
        
        
        //二次项 (M_i*x_i1)/(2*h_i)-(M_i1*x_i)/(2*h_i)
        //     =(M_i*x_i1-M_i1*x_i)/(2*h_i)
        
        arb_mul(s,arb_mat_entry(M,j,0),array_x[i+1],prec);
        arb_mul(w,arb_mat_entry(M,j+1,0),array_x[i],prec);
        arb_sub(s,s,w,prec),
        
        arb_mul_ui(w,h_i,2,prec); //最后相除
        arb_div(arry_fit_save[i].b,s,w,prec);
        
        //三次项 M_i1/(6*h_i)-M_i/(6*h_i)
        //     =(M_i1-M_i)/(6*h_i)
        
        arb_sub(s,arb_mat_entry(M,j+1,0),arb_mat_entry(M,j,0),prec);
        
        arb_mul_ui(w,h_i,6,prec); //最后相除
        arb_div(arry_fit_save[i].a,s,w,prec);
        
    }
    
    arb_clear(s); //释放变量
    arb_clear(t);
    arb_clear(w);
    arb_clear(h_i);
    
    printf("多项式系数求解完成\n\n");
    
    /*
    //将三次样条拟合函数的整个结构以二进行形式存下来
    FILE *fitp ;
    fitp = fopen(fitted_file,"wb"); // b:表示以二进制写入
    
    fwrite(arry_fit_save, sizeof(struct FUNC_FITTED_DATE), entries_number, fitp); //2：表示将数组中两个元素写入文件
    
    fclose(fitp);
    */
    
    //暂时不用二进制的形式保存数据，因为保存的数据不能单独读取
    //猜测与 arb_t 的特殊结构有关 
    //采用文本文件的形式
    FILE* fitp;
    fitp = fopen(fitted_file, "w");
    
    if (!fitp) {
        printf("Math_Method -> fit_date -> ERROR: Unable to open file: %s", fitted_file);
        exit(1);
        return EXIT_FAILURE;
    }
    
    for(slong i=entries_number-1; i >= 0; i=i-1 ) //倒序输出，读入还会倒一下，就正了
    {
        //存储， 一个点 (a,c,d,x) 一行，用 tab 分隔
        
        arb_dump_file(fitp,arry_fit_save[i].a); 
        fprintf(fitp, "\t");
        arb_dump_file(fitp,arry_fit_save[i].b);
        fprintf(fitp, "\t");
        arb_dump_file(fitp,arry_fit_save[i].c); 
        fprintf(fitp, "\t");
        arb_dump_file(fitp,arry_fit_save[i].d);
        fprintf(fitp, "\t");
        arb_dump_file(fitp,arry_fit_save[i].x);
        fprintf(fitp, "\n");
    }
    fclose(fitp);
    
    printf("已拟合结果将存到： %s\n\n", fitted_file);
    
    //全部完成，释放变量
    free(array_x);
    free(array_y);
    free(arry_fit_save);
    arb_mat_clear(A);
    arb_mat_clear(B);
    arb_mat_clear(M);
    
    //printf("完成数据存储\n\n");
    
    return 0;
}

//拟合数据拟合恢复
int Func_fit_restore(struct FUNC_FITTED_DATE **fit_res_2, char* fitted_file, slong prec)
{
    /*
    //二进制读取，暂时不可用
    ulong my_number;
    struct stat temp_buf; //利用 stat 获取相应数组的大小， stat 可以获取文件的各种各样的信息
    
    stat(fitted_file, &temp_buf);
    my_number = temp_buf.st_size/sizeof(struct FUNC_FITTED_DATE);
    
    //free(temp_buf);
    //根据所获得的点数，用动分配内存
    fit_res=(struct FUNC_FITTED_DATE *) calloc(my_number,sizeof(struct FUNC_FITTED_DATE));
    
    //arb_t 使用时，需初始化，即使是数组，也得初始化
    for(ulong i=0; i < my_number; i+=1)
    {
        arb_init(fit_res[i].a);
        arb_init(fit_res[i].b);
        arb_init(fit_res[i].c);
        arb_init(fit_res[i].d);
        arb_init(fit_res[i].x);
    }
    
    
    FILE *fp ;
    
    fp = fopen(Func_output_fitted_file,"rb");
    
    fread(fit_res,sizeof(struct FUNC_FITTED_DATE),my_number,fp);
    
    fclose(fp);
    */
    
    
    //文本文件读取
    //定义一个由 FUNC_FITTED_DATE 构成的单链列表
    SLIST_HEAD(func_fitted_list, FUNC_FITTED_DATE) head = SLIST_HEAD_INITIALIZER(head);
    SLIST_INIT(&head);
    
    printf("#### 开始读取拟合结果 ####\n");
    
    //读取还应的数据文件
    FILE* input_file = fopen(fitted_file, "r"); //打开文件
    
    if (!input_file) {
        printf("ERROR: Unable to open file: %s", fitted_file);
        
        return EXIT_FAILURE;
    }
    
    ulong row = 0; //有效行数，计算成功读取的行数，里面可能有空行等情况
    ulong row_s = 0; //总行数，用以把最开始的几行注释忽略
    
    char* str_a="";
    char* str_b="";
    char* str_c="";
    char* str_d="";
    char* str_x="";
    
    //然后逐行读取文件，直到出现错误或文件结束
    while (!ferror(input_file) && !feof(input_file))
    {
        size_t buffer_size = 0; //size_t是无符号的长整型
        //一般用来表示字节数的多少。常用于如sizeof返回值的类型
        //size_t是一种跟具体的平台有关联的类型，会具体调整其能表示的范围，因此其可移植性会更好
        //size_t能保证可以存储任何类型理论上可能的对象的最大值，包括数组类型。
        char *buffer = NULL;
        
        //getline() 函数可读取文件中的整行，并负责分配必要的内存
        getline(&buffer, &buffer_size, input_file);
        
        
        //读取整行后，使用 strtok() 函数将每一行分成字元「token」。遍历字元，选择所需的列
        //这里理应 strlen(buffer) > 0 但最后一行会遇到一个奇怪的字符，导至多出一行
        //现暴力的将其设为 > 6 个字符长度 
        // 空行 buffer[0] !='\n' ，用前面6个字符长度已不需要
        if (row_s >= Func_output_skip_header && strlen(buffer) > 6)
        {
            unsigned int column = 0; //对于列数来说，足够
            
            char* token = strtok(buffer, Func_output_delimiter);
            
            while (token != NULL) //读取一行中的各例，读取的内容是字符串的形式
                //可通过 sscanf(token, "%lf", &value); 的形式
                //把字符串的内容给 double value 变量
            {
                //这里输入的是字符串
                //由 int arb_load_str(arb_t x, const char *str) 直接读入即可
                if (column == 0) {
                    str_a = token;
                }else if (column == 1) {
                    str_b = token;
                }else if (column == 2) {
                    str_c = token;
                }else if (column == 3) {
                    str_d = token;
                }else if (column == 4) {
                    str_x = token;
                }
                
                column += 1;
                token = strtok(NULL, Func_output_delimiter);
            }
            
            //当选择了 x 和 y 值时，将新数据点插入链表中：
            //malloc() 函数为新数据点动态分配（保留）一些持久性内存
            //struct FUNC_FITTED_DATE *datum = malloc(sizeof(struct FUNC_FITTED_DATE));
            struct FUNC_FITTED_DATE *datum = (struct FUNC_FITTED_DATE *)calloc(1,sizeof(struct FUNC_FITTED_DATE));
            
            int jude_1,jude_2,jude_3,jude_4,jude_5;
            
            //arb_t 使用时，需初始化，即使是结构体，也得初始化
            arb_init(datum->a); //arb_t 使用前初始化，避免可能存在的内存分配问题
                                //如没初始化时，当其遇是整数的时候，不会按 arb_t 结构体 48 分配空间
                                //这样，就导制存存空间分配的问题，发生段错误 Segmentation fault (core dumped)
            arb_init(datum->b);
            arb_init(datum->c);
            arb_init(datum->d);
            
            jude_1=arb_load_str(datum->a,str_a);
            jude_2=arb_load_str(datum->b,str_b);
            jude_3=arb_load_str(datum->c,str_c);
            jude_4=arb_load_str(datum->d,str_d);
            jude_5=arb_load_str(datum->x,str_x);
            
            if( jude_1 != 0 || jude_2 != 0 || jude_3 != 0 || jude_4 != 0 || jude_5 != 0)
            {
                printf("Math_Method -> fit_date -> 加入数据失败，行数：%li，文件名：%s\n",row_s,Func_output_file);
                exit(1);
            }
            
            SLIST_INSERT_HEAD(&head, datum, func_fitted_date_entries);
            
            row+= 1; //有效行数
        }
        
        //释放内存，改变行数，读取下一行
        free(buffer);
        row_s += 1; //总行数
    }
    
    fclose(input_file); //读取文件完成，关闭文件
    
    if(row==0)
    {
        printf("读取错误，条目为 0 ：%s\n",fitted_file);
        return EXIT_FAILURE;
    }
    
    
    //由于你将不知道要创建的数组的大小，因此必须手动分配它们的内存
    //注意，这里分配的内存，即数组中变量的个数，需与数据的个数严格匹配
    
    const size_t entries_number = row; //计算所有数据行数，得到娈量个数
    
    //将相应的数据用存到数组中，供下一步使用
    //手动分配内存的方法有两种 calloc 和 malloc ，malloc在数据点很少「如10」时有问题
    struct FUNC_FITTED_DATE * fit_res = (struct FUNC_FITTED_DATE *) calloc(entries_number,sizeof(struct FUNC_FITTED_DATE)); //根据变量的数目分配内存
    
    //fit_res =(struct FUNC_FITTED_DATE *) malloc(entries_number * sizeof(struct FUNC_FITTED_DATE)); //根据变量的数目分配内存
    
    
    //arb_t 使用时，需初始化，即使是数组，也得初始化
    for(ulong i=0; i < entries_number; i+=1)
    {
        arb_init(fit_res[i].a);
        arb_init(fit_res[i].d);
        arb_init(fit_res[i].c);
        arb_init(fit_res[i].d);
        arb_init(fit_res[i].x);
    }
    
    
    
    if ( !fit_res ) {
        printf("ERROR: Unable to allocate data arrays\n");
        return EXIT_FAILURE;
    }
    
    
    struct FUNC_FITTED_DATE *datum;
    datum = SLIST_FIRST(&head);
    
    ulong i = 0; //给数组 x y 赋值用
    
    //遍历链表，存入数组
    SLIST_FOREACH(datum, &head, func_fitted_date_entries)
    {
        arb_set(fit_res[i].a,datum->a); //存值到数组
        arb_set(fit_res[i].b,datum->b);
        arb_set(fit_res[i].c,datum->c);
        arb_set(fit_res[i].d,datum->d);
        arb_set(fit_res[i].x,datum->x);
        i += 1;
    }
    
    //处理完链表后，需清理
    //手动释放已分配的内存，防止内存泄漏
    while (!SLIST_EMPTY(&head))
    {
        struct FUNC_FITTED_DATE *datum = SLIST_FIRST(&head);
        
        SLIST_REMOVE_HEAD(&head, func_fitted_date_entries);
        
        free(datum);
    }
    Fit_entries_number=entries_number; //向全局变量传递整个数据点的多少
    
    
    *fit_res_2=fit_res;  //把二级指针的内容修改为 fit_res 所指向的地址
                        //而这个二级指针指向的是一个外部指针的地址
                        //从而将外部指针指向的地址修改
    
    printf("\n读取成功： %s\n有效数(a,b,c,d,x) 共计：%li\n\n", fitted_file,entries_number);
    
    
    
    //printf("读取完成，共计点数： %li\n", entries_number);
    
    return 0;
    
}


//求相应拟合函数的值
int Fit_get_value(arb_t res, const arb_t x, const struct FUNC_FITTED_DATE* Fit_func,
                  const int order, slong prec)
{
    arb_t s;
    arb_init(s);
    
    //判断x所属区间
    if ( arb_lt(x,Fit_func[0].x) || arb_gt(x,Fit_func[Fit_entries_number-1].x) )
    {
        printf("Math_Method -> fit_date -> 所给参数不在拟合区间内\n");
        exit(1);
    }
    
    
    //对x查找所属区间
    //算法后面优化
    ulong my_i=-1;
    for (ulong i=0; i <  Fit_entries_number-1; i+=1)
    {
        if (arb_lt(x,Fit_func[i+1].x) )
        {
            //printf("x: ");arb_printd(x,7);printf("\n");
            //printf("res: ");arb_printd(Fit_func[i].x,7);
            //printf("\t");arb_printd(Fit_func[i+1].x,7);printf("\n");
            
            my_i=i;
            
            break;
        }
    }
    
    //输出拟合函数值
    if( order==0 ) //原函数
    {
        //a*x^3 + b*x^2 + c*x + d
        arb_sqr(res,x,prec); //三次
        arb_mul(res,res,x,prec);
        arb_mul(res,res,Fit_func[my_i].a,prec);
        
        arb_sqr(s,x,prec); //二次
        arb_mul(s,s,Fit_func[my_i].b,prec);
        arb_add(res,res,s,prec);
        
        arb_mul(s,x,Fit_func[my_i].c,prec); //一次
        arb_add(res,res,s,prec);
        
        arb_add(res,res,Fit_func[my_i].d , prec);
        
    }else if ( order==1 ) //一阶函数
    {
        //3*a*x^2 + 2*b*x + c 
        
        arb_sqr(s,x,prec); //二次
        arb_mul(s,s,Fit_func[my_i].a,prec);
        arb_mul_si(res,s,3,prec);
        
        arb_mul(s,x,Fit_func[my_i].b,prec); //一次
        arb_mul_ui(s,s,2,prec);
        arb_add(res,res,s,prec);
        
        arb_add(res,res,Fit_func[my_i].c, prec);
        
    }else if ( order==2 ) //二阶函数
    {
        //6*a*x + 2*b
        
        arb_mul(s,x,Fit_func[my_i].a,prec); //一次
        arb_mul_ui(res,s,6,prec);
        
        arb_mul_si(s,Fit_func[my_i].b,2,prec);
        arb_add(res,res,s, prec);
    }
    
    
    arb_clear(s);
    
    return 0;
}

