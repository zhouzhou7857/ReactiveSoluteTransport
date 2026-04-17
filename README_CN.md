# ReactiveSoluteTransport

## 项目概述

ReactiveSoluteTransport 是一个基于离散裂隙网络（DFN）的二维溶质运移模拟程序，主程序位于 `Code/src`。当前主线流程是：

1. 读取域、模拟、DFN 和 chemistry 输入文件；
2. 构建或生成 DFN；
3. 计算裂隙网络流动并提取 backbone；
4. 执行粒子追踪；
5. 按粒子携带反应物的损耗更新裂隙 aperture；
6. 在 aperture 变化后重建网络并重算流场；
7. 输出突破曲线、粒子快照和 DFN 演化结果。

## 当前主线特点

- 以粒子追踪方式模拟裂隙内平流与基质扩散。
- 支持无限基质与有限基质两种输运模式。
- 支持均匀注入、局部注入和入口通量加权注入。
- 当前主路径中的 aperture 更新是 chemistry-driven，而不是旧的粒子通过量经验公式。
- 支持从 `Input/Chemistry_files/` 读取化学参数。
- 支持 `generation_realistic3` 和 `generation_realistic4` 的构造函数版 DFN 生成路径。
- 运行时会额外输出 `Output/DFN_raw.txt`，用于保存 backbone 之前的原始 DFN。

## 技术栈

- 语言：C++
- 几何：CGAL
- 数值与特殊函数：Boost
- 线性代数/求解：uBLAS 风格结构与 SuiteSparse 相关依赖
- 可视化：OpenGL / GLUT
- 工程方式：Eclipse CDT + GNU Make

## 主流程

### 输入读取

`Input/File_names.txt` 当前按 4 行读取：

1. 域文件
2. 模拟文件
3. DFN 文件
4. chemistry 文件

`Parameters` 会进一步读取：

- `Input/Domain_files/` 中的域尺寸、基质扩散系数、孔隙度和左右边界水头；
- `Input/Simulation_files/` 中的粒子数、时间范围、反应时间步、输出时间步等；
- `Input/DFN_files/` 中的 DFN 模式与参数；
- `Input/Chemistry_files/` 中的初始反应物浓度、衰减率、化学计量系数、矿物摩尔体积和裂隙厚度。

### DFN 生成与读取

当前主程序在 `Code/src/PERFORM.cpp` 中分两类处理 DFN：

- `generation_realistic3` / `generation_realistic4`
  - 走构造函数版：
    `NetworkMeshes(param.density_param, param.exponent_param, param, domain)`
- 其他模式
  - 仍走原来的文件/模式构造函数：
    `NetworkMeshes(param.code_path, param.file_name_DFN, domain)`

对于 `generation_realistic3`，当前主路径下：

- 以 `density_param` 作为累计密度停止条件，而不是精确裂缝数；
- 用 `exponent_param` 控制长度分布；
- 用 `b_min`、`b_max`、`mean_lnb`、`RSD_lnb` 控制 aperture 分布；
- 当前实现中长度与 aperture 采用同一个随机数驱动，因此是正相关设置。

### Chemistry-driven aperture 更新

当前主路径下，几何更新采用以下链条：

1. 每个新注入粒子被赋予 `INITIAL_REACTIVE_CONCENTRATION`；
2. 粒子在裂隙段中的停留时间导致浓度衰减；
3. 浓度损失乘以粒子代表体积得到反应物摩尔损失；
4. 再经化学计量关系和矿物摩尔体积换算为矿物体积变化；
5. 矿物体积变化换算为 segment 平均 aperture 增量；
6. aperture 更新后，若有裂隙关闭，则重建网络并重算流场。

旧的 `ChangeAperture()`、`ComputeAperture()` 等经验接口仍保留在代码里，但当前主路径默认不使用。

## 核心文件

- `Code/src/PERFORM.cpp`
  - 主入口，负责参数读取、DFN 构建分流、backbone、输运和结果输出。
- `Code/src/Input_Output/Parameters.cpp`
  - 读取域、模拟、DFN 和 chemistry 参数。
- `Code/src/Domain_Definition/NetworkMeshes.cpp`
  - DFN 读取、DFN 生成、交点处理、backbone、流速计算。
- `Code/src/Transport/Transport.cpp`
  - 主输运循环、粒子注入、步进输运、chemistry 耦合 aperture 更新。
- `Code/src/Chemistry/Chemistry.cpp`
  - 反应物浓度衰减与 aperture 增量换算。
- `Code/src/Input_Output/Results.cpp`
  - 输出 `cdf.txt` 和 `pdf.txt`。

## 当前输入示例

仓库中当前提供了以下与 `generation_realistic3` 相关的示例文件：

- `Input/Domain_files/Domain_10cm_square.txt`
- `Input/DFN_files/DFN_realistic3_10cm_target200.txt`
- `Input/DFN_files/DFN_generator_description.txt`

这些文件用于说明 10 cm × 10 cm 方形域和 `generation_realistic3` 参数输入方式。

## 输出文件

当前主流程运行后常见输出包括：

- `Output/DFN_raw.txt`
  - backbone 之前的原始 DFN
- `Output/DFN_init.txt`
  - backbone 后、输运更新前的初始 DFN
- `Output/DFN.txt`
  - 当前/最终 DFN
- `Output/DFN_aperture_delta.txt`
  - 相对初始 DFN 的 aperture 变化
- `Output/DFN_step*.txt`
  - 随反应时间步记录的 DFN 快照
- `Output/particle_positions_t*.csv`
  - 粒子位置快照
- `Output/cdf.txt`
  - 到达时间累计分布
- `Output/pdf.txt`
  - 到达时间概率密度

## 当前代码说明

- 当前主线说明以 `Code/src` 为准。
- `Code/Backup0318`、`Code/Backup0318_unmodified` 仅为历史参考。
- `Code/Release` 和 `Code/Debug` 是构建产物目录，不应作为主逻辑说明来源。
