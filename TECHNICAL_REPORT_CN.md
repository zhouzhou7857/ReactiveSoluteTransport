# 基于实现的技术报告

## 说明

本报告基于当前有效代码 `Code/src` 与当前输入格式 `Input/` 编写。

## A. 总体目标

该代码用于模拟二维离散裂隙网络（DFN）中的粒子追踪溶质运移，并将输运过程与 DFN 结构演化耦合。当前主路径中的耦合方式是 chemistry-driven aperture 更新：

- 每个注入粒子携带初始反应物浓度；
- 粒子在裂隙段中停留时，反应物浓度随时间衰减；
- 反应物损耗换算为矿物体积变化；
- 矿物体积变化换算为裂隙段平均 aperture 增量；
- aperture 变化后，裂隙网络流场按需要重新计算。

这与旧版本基于粒子通过量的经验 aperture 更新不同。旧接口仍保留在代码中，但不再是当前主执行路径的默认机制。

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
2. `Parameters::Parameters()` 读取 `Input/File_names.txt`
3. `Parameters::read_param()` 读取：
   - 域参数
   - 模拟参数
   - DFN 模式与 DFN 参数
   - 可选 chemistry 参数
4. `PERFORM.cpp` 通过 `ConfigureChemistryParameters(...)` 配置当前 chemistry 常量
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

## C. DFN 构建路径

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

## D. 粒子模型

每个粒子保存以下关键信息：

- 当前位置
- 当前网格编号
- 前一个网格编号
- 网格与交点访问历史
- 当前时间
- 注入时间
- 当前裂隙内累计平流距离
- 当前裂隙内累计停留时间
- 粒子当前携带的反应物浓度
- 粒子代表的流体体积

粒子注入时会被赋予：

- 入口位置和入口裂隙
- 清零后的局部输运状态
- `reactive_concentration = INITIAL_REACTIVE_CONCENTRATION`
- 按入口通量和粒子时间份额计算得到的 `representative_volume`

## E. Chemistry-driven 几何更新

### 当前主路径的更新方式

当前主路径中的 aperture 更新不再依赖旧的 `Nb/Nt` 经验公式。

当前链条是：

1. `UpdateReactiveConcentration()` 计算：
   - `C_out = C_in * exp(-k * t_residence)`
2. `EvaluateReactiveStep()` 继续计算：
   - 反应物摩尔损失
   - 矿物摩尔变化
   - 矿物体积变化
   - 局部 aperture 变化
   - 段平均 aperture 变化
3. `Transport.cpp` 将这些变化累计到 `step_aperture_change`
4. 当前裂隙段 aperture 加上该段平均增量
5. 若 aperture 归零，则可能触发网络重建和流场重算

### 所需 chemistry 参数

chemistry 文件当前提供：

- 初始反应物浓度
- 浓度衰减率
- 反应物到矿物的化学计量系数
- 矿物摩尔体积
- 裂隙面外厚度

## F. 旧接口的状态

以下旧接口仍保留在代码中：

- `ComputeAperture()`
- `ComputeDeltaAperture()`
- `NetworkMeshes::ChangeAperture()`
- `NetworkMeshes::ChangeApertureByRatio()`

它们主要用于兼容性和历史参考，不是当前主执行路径的默认机制。

## G. 输出

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
