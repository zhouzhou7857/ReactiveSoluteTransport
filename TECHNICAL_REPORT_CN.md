# 基于实现的技术报告

## 说明

本报告基于当前有效代码 `Code/src` 与当前维护的输入格式 `Input/` 编写。

## A. 总体目标

该代码用于模拟二维离散裂隙网络（DFN）中的粒子追踪溶质运移，并将输运过程与裂隙几何演化耦合。当前主路径中的耦合方式已经变成基于 PHREEQC 标定的累计矿物体积变化律：

- 每个注入粒子携带一个代表性流体体积；
- 每个粒子还保存累计反应年龄 `particle_age`；
- chemistry 模块计算参考累计律 `DeltaV_ref(t)`；
- 单个输运时间段上的矿物体积增量为
  `DeltaV_ref(t_end) - DeltaV_ref(t_start)`；
- 该体积增量再乘以 `V_particle / Vref`，并可选乘以额外几何缩放因子；
- 各 segment 上累计的矿物体积变化再换算为平均 aperture 增量；
- aperture 更新后按需要重算流场。

旧的浓度衰减主路径和 `Nb/Nt` 经验 aperture 公式仍保留在代码中，但不再是当前默认执行路径。

## B. 主流程

### 入口与主要模块

主入口：

- `Code/src/PERFORM.cpp`

主要模块：

- `Code/src/Input_Output/Parameters.h/.cpp`
- `Code/src/Domain_Definition/NetworkMeshes.h/.cpp`
- `Code/src/Domain_Definition/DFNComputation.h/.cpp`
- `Code/src/Transport/Transport.h/.cpp`
- `Code/src/Chemistry/Chemistry.h/.cpp`
- `Code/src/Input_Output/Results.h/.cpp`

### 当前执行顺序

1. `PERFORM.cpp` 创建 `Parameters param`
2. `Parameters::Parameters()` 读取 `Input/File_names.txt`，或读取命令行传入的替代文件列表
3. `Parameters::read_param()` 读取：
   - 域参数
   - 模拟参数
   - DFN 模式与 DFN 参数
4. `PERFORM.cpp` 先重置 chemistry 默认值，再根据：
   - simulation input 中的 `Vref` 和 `fracture thickness`
   - 运行时环境变量 `RST_CHEM_MODE`、`RST_DELTA_V_A*`、`RST_DELTA_V_K*`、`RST_DELTA_V_L`、`RST_USE_EFFECTIVE_DIFFUSION_HEIGHT_FACTOR`、`RST_USE_VP_WIDTH_CORRECTION`
   配置当前活跃 chemistry
5. `PERFORM.cpp` 根据 DFN 模式分流：
   - `generation_realistic3/4`：
     `NetworkMeshes(param.density_param, param.exponent_param, param, domain)`
   - 其他模式：
     `NetworkMeshes(param.code_path, param.file_name_DFN, domain)`
6. 在 backbone 之前输出 `Output/DFN_raw.txt`
7. 若 `backbone=true`，调用 `return_backbone()` 提取连通流动骨架
8. `Transport::Particles_Transport()` 执行步进输运与几何更新
9. `Results` 输出 `cdf.txt` 与 `pdf.txt`
10. 程序写出：
   - `DFN_init.txt`
   - `DFN.txt`
   - `DFN_aperture_delta.txt`

## C. 输入格式

### File list

`Input/File_names.txt` 现在主用为 3 行格式：

1. domain file
2. simulation file
3. DFN file

第 4 行 chemistry file 仍然兼容旧格式，但当前主路径不会再读取 chemistry input 文件。

### Simulation 文件

simulation 文件末尾现在包含两个真正参与当前 chemistry 主路径的参数：

11. `Vref` `[m^3]`
12. `fracture thickness` `[m]`

它们在 `Parameters.cpp` 中读入，并在 `PERFORM.cpp` 中传给 chemistry 模块。

## D. DFN 构建路径

### 从旧版 `realistic3` 到当前版本的切换

这次切换的核心机制有两步：

1. 在 `Parameters.cpp` 中，把 `generation_realistic3/4` 需要的参数真正读入 `Parameters`
2. 在 `PERFORM.cpp` 中，对 `generation_realistic3/4` 分流到构造函数版
   `NetworkMeshes(double,double,Parameters,Domain)`

也就是说，当前主程序不再让 `generation_realistic3` 默认走旧的文件驱动
`NetworkMeshes(code_path,file_name_DFN,domain)` 路径，而是切到了构造函数版生成器。

### 当前 `generation_realistic3` 的行为

当前有效实现位于 `NetworkMeshes.cpp` 的构造函数分支中。

其特点是：

- 以 `density_param` 作为累计密度停止条件
- 裂隙中心在域内均匀随机
- 裂隙长度由 `exponent_param` 控制
- 裂隙方向在 `[0, pi]` 上均匀随机
- aperture 来自截断对数正态分布，参数为：
  - `b_min`
  - `b_max`
  - `mean_lnb`
  - `RSD_lnb`
- 当前实现中长度和 aperture 使用同一个随机数驱动，因此更长的裂隙倾向于更大的 aperture

注意：

- `density_param` 不是精确裂缝条数
- 当前 `generation_realistic3` 是“密度阈值停止”，不是“固定裂缝条数停止”

## E. 粒子模型

每个粒子保存以下关键信息：

- 当前位置
- 当前网格编号
- 前一个网格编号
- 网格与交点访问历史
- 当前时间
- 注入时间
- 当前裂隙内累计平流距离
- 当前裂隙内累计停留时间
- 累计 chemistry 年龄 `particle_age`
- 粒子代表的流体体积 `representative_volume`

粒子注入时会被赋予：

- 入口位置和入口裂隙
- 清零后的局部输运状态
- 按入口通量和粒子时间份额计算得到的 `representative_volume`：
  `V_particle = |u_inlet| * aperture_inlet * thickness * dt_particle`

## F. 当前活跃 chemistry-driven 几何更新

### 当前主路径的更新方式

当前主路径中的 aperture 更新基于 `Chemistry.cpp` 中的累计矿物体积变化律：

`DeltaV_ref(t) = A1 * (1 - exp(-k1 * t)) + A2 * (1 - exp(-k2 * t)) + L * t`

单个粒子在一个输运时间段上的体积增量为：

`DeltaV_particle[t_start,t_end] = (DeltaV_ref(t_end) - DeltaV_ref(t_start)) * scale`

其中 `scale`：

- 至少包含 `V_particle / Vref`
- 可选乘以 effective diffusion-height factor
- 可选乘以 Vp-width correction

`Transport.cpp` 会先把每个 segment 上的矿物体积变化累计起来，再统一执行：

`delta_b = delta_V / (segment_length * thickness * 2)`

### 可选缩放控制

当前代码支持以下可选运行时控制：

- `RST_USE_EFFECTIVE_DIFFUSION_HEIGHT_FACTOR`
- `RST_EFFECTIVE_DIFFUSION_COEFFICIENT`
- `RST_EFFECTIVE_DIFFUSION_TIME`
- `RST_USE_VP_WIDTH_CORRECTION`

这些主要用于 benchmark 和敏感性分析，不是默认必须开启的项。

## G. 旧接口状态

以下旧接口仍保留在代码中：

- `UpdateReactiveConcentration()`
- `EvaluateReactiveStep()`
- `ComputeAperture()`
- `ComputeDeltaAperture()`
- `NetworkMeshes::ChangeAperture()`
- `NetworkMeshes::ChangeApertureByRatio()`

它们主要用于兼容性和历史参考，不是当前主执行路径的默认机制。

## H. 输出

当前运行可能生成：

- `Output/DFN_raw.txt`
- `Output/DFN_init.txt`
- `Output/DFN.txt`
- `Output/DFN_aperture_delta.txt`
- `Output/DFN_step*.txt`
- `Output/particle_positions_t*.csv`
- `Output/cdf.txt`
- `Output/pdf.txt`

其中 `DFN_raw.txt` 用于查看 backbone 筛选前的原始 DFN。
