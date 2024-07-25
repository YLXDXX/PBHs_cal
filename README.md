# 程序说明

## 简介

这是一个原初黑洞（PBHs）丰度及标量诱导引力波计算程序，为保证计算准确性，采用了 FLINT 中的任意精度计算模块 `arb`，为保证计算速度，采用 C 语言实现

有两种方法来计算PBHs的丰度：

1. 基于 compaction function 的 Press–Schechter 方法
2. 基于 compaction function 的 Peak Theory 的方法

有三种方法计算高斯情况下的标量诱导引力波，有一种方法计算非高斯情况下的标量诱导引力波「 $f_{NL}$ 」



### 框架

其中 `main.c`​ 文件为主程序入口，可在其中可进行各种量的设定，计算相关参数设定在 `src/set` ​下，具体请查看源码及相关注释

而 `main.c`​ 中的相关计算，又是调用如下五个模块来完成的「在 `src/func`​ 下」：

1. General

    包含一些基础通用的计算和设定，供后面 PS 和 PK 调用
    
2. PS_Method

    PS 方法的计算相关
    
3. Peak_Theory

    Peak Theory 计算相关
    
4. GW_induced

    标量诱导引力波计算相关

4. Math_Method

    用 `arb`​​ 实现所需的各种数学算法，供前面各计算调用
    
5. Other

    用于一些程序中其它功能的实现
    
6. Test

    包含一些测试相关的内容，用于验证各计算模块是否正确

另外：

`new_type.h`​ 定义了 PBHs 计算所需的新类型

`phy_constant.c` 和 `phy_constant.h` 定义了 PBHs 计算所需的各种变量，这里的量都是作为全局变量的形式存在，一个全局变量需要在这两个文件中都进行设定，如果是类似于 `arb_t` 这样需要初始化的全局变量，统一在 `src/set/00_global_variable.c` 中进行初始化。

> 注意，六个大模块中各小模块的顺序是按该模块的实际计算顺序排列
>
> 而 `src/set` ​中各模块的设定顺序，也是按照实计算顺序需要排列

‍



## 计算模块介绍

### General

具体包含：

`basis.h`​ 最基本的初始设定：功率谱，窗口函数，转移函数

`non_gaussianity.h`​ 局域非高斯性 $\zeta=f(\zeta_G)$

`typical_profile.h`​ 计算 typical profile $\hat{\zeta}(r)$

`compaction_func.h`​ 计算 compaction function $\mathcal{C}(r)$ 及其第一个极大值

`threshold.h`​ 计算阈值 $\mathcal{C}_{th}$ 和 $\mathcal{C}_{\ell,th}$

‍`other_func.h` 在PS和PT方法中需要的一些其它计算



### Peak Theory

具体包含：

`generate_mass.h`​ PBHs 生成质量与当时视界质量间的关系

`number_density.h`​ 计算峰的数密度

`abundance.h`​ 计算 PBHs 丰度

‍

### PS Method

具体包含如下：

`covariance.h`​ 计算协方差 $\Sigma_{XX},\Sigma_{XY},\Sigma_{YX},\Sigma_{YY}$

`zeta_probability.h`​ 计算 $P(\zeta)$，包含非高斯性概率的计算

`C_l_probability.h`​​ 计算线性的 compaction function 的概率 $P(\mathcal{C}_\ell)$

`ps_generate_mass.h`​ PBHs 生成质量与当时视界质量间的关系

`ps_abundance.h`​ 计算 PBHs 丰度

`ps_abundance_all_k.h`​ 计算 PBHs 丰度，主要针对连续谱，考虑了连续谱中所有 $k$ 模式的影响



### GW_induced

具体包含如下：

`**_help_func.h` 计算诱导引力波辅助函数

`**_power_spectra.h`  计算诱导引力波的功率谱和能量密度谱

`gw_induced.h`  通用计算诱导引力波的功率谱和能量密度谱函数（高斯情况）

`03_power_spectra.h` 通用计算诱导引力波的功率谱和能量密度谱函数（非高斯情况）



### Math Method

具体包含如下：

`func_constant.h`​​ 一些需要储存的数据，如 gauss_kronrod 积分的系数和权重

`find_root.h`​ 求根算法，包含 Brent 方法和折半查找法，可查找某区间内的单个根和多个根

`fit_date.h`​ 数据拟合相关，用于化简 PK 相关计算

`quadrature.h`​ 一元积分算法，可采用：自适应 Simpson 积分、自适应 Gauss–Kronrod 积分、Double Exponential 积分，其中 Gauss–Kronrod 积分有迭代和递归两个版本

`quadrature_binary.h ` 二元定积分算法，包括矩形边界条件和由函数描述的边界条件，可采用 Gauss–Kronrod 积分或 Double Exponential 积分。对于二维矩形积分区域，有自适应和非自适应两个版本，对于 Gauss–Kronrod 积分而言，自适应比非自适应快，对于 Double Exponential 而言，非自适应远快于自适应，为提高精度，可对非自适应提高迭代次数。

另外，对于其中的积分计算，引入了多线程方式，但需要注意的是，并不是线性量越多越快，需要根据相对积分函数的性质灵活选取。



### Other

具体包含如下：

`fitting_DOF.h` 获取不同温度 $T$ 或不同波数 $k$ 对应的自由度数

`output.h`​ 用于获取存储路径、输出数据的进度条展示等

`other.h` 诱导引力波频率 $f$ 与波数 $k$ 的转换等



### Test

具体包含如下：

`func_test.h`​ 定义各种测试相关函数

`routine_test.h`​ 相关测试的具体实现





## 相关参数设定

计算相关参数，设定文件为 `main.c`​ 和 `src/set`​ 下的文件，其中 `main.c`​ 是一些宏观量的设定，而 `src/set` ​的设定较为具体：

`00_global_variable.c`​ 各个全局变量的初始化「有一些极其特殊函数中的全局变量，不需要在这里初始化」

`01_cal_math.c`​ 常数和数学计算行为

* 对于 broken power law 类型的功率谱，建议选择 `gauss_kronrod_iterate`​ 积分方式
* 对于 log-normal 类型功率谱，建议选择 `double_exponential`​ 积分方式

`02_physical_parameter.c`​ 宇宙学中的基本参数

`03_power_spectra.c`​ 功率谱的具体参数

`04_zeta_r_cal.c`​ 利用 mean profile 计算 $\zeta(r)$

`05_main_cal.c`​ 本计算程序的主要计算参数设定，包括：非高斯性参数，各种功率谱和各种非高斯性情况下 $r_m$，$\mu_{th}$，$\beta$，$f$ 等计算参数，临界坍缩、诱导引力波计算等

`06_draw_pic.c`​ 结果输出，用于画图





## 项目编译

本项目采用 [xmake](https://xmake.io/) 进行构建，依赖：

- [FLINT](https://flintlib.org/) 中的 `arb`​ 模块，建议安装使用 FLINT 3.x series
- [Cuba](https://feynarts.de/cuba/) 中的多维积分模块 <cubaq.h>

> 注意，FLINT 2.x series 不包含 `arb`​ 模块

> 注意，在编译 Cuba 时，需加上参 --with-real=16，这样才可以用 <cubaq.h> 模块，不然是用的 <cuba.h> 模块



安装好 Flint 后，需修改 `xmake.lua`​ 文件，将中的 `add_includedirs`​ 选项改为自己 Flint 安装的地址。然后在 `xmake.lua`​ 所在文件夹下进行编译

```shell
xmake b
```

编译完成后运行程序

```shell
xmake r
```

使用环境示例：

* 系统：Debian 12
* FLINT： 3.1.2
* Cuba:：4.2.2
* 编译工具： Gcc 12.2.0
* 构建工具： Xmake 2.8

‍

### 注：

- `wxMaxima`​文件夹内的东西，是在编写程序时借助 wxMaxima 所做的一些数学运算

‍
