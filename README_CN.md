# ReactiveSoluteTransport

## 项目概述

ReactiveSoluteTransport 是一个基于离散裂隙网络（DFN）的二维溶质运移模拟程序，主程序位于 `Code/src`。当前主线流程是：

1. 读取域、模拟和 DFN 输入文件；
2. 构建或生成 DFN；
3. 计算裂隙网络流动并提取 backbone；
4. 执行粒子追踪；
5. 按粒子年龄驱动的 `DeltaV` 体积增量更新裂隙 aperture；
6. 在 aperture 变化后重建网络并重算流场；
7. 输出突破曲线、粒子快照和 DFN 演化结果。

## 当前主线特点

- 以粒子追踪方式模拟裂隙内平流与基质扩散。
- 支持无限基质与有限基质两种输运模式。
- 支持均匀注入、局部注入和入口通量加权注入。
- 当前主路径中的 aperture 更新是基于 `DeltaV` 累积律的 chemistry-driven 路径，而不是旧的粒子通过量经验公式。
- 当前 V5-final benchmark 主要面向过饱和水中的 gypsum precipitation，以及类似的“已知参考反应曲线 + 粒子代表体积权重”的沉淀/溶解问题。
- 当前主路径不再从 `Input/Chemistry_files/` 读取活跃 chemistry 参数。
- `Vref` 和 fracture thickness 由 simulation input 提供，其余 `DeltaV` 系数默认来自代码或运行时环境变量覆盖。
- 支持运行时 chemistry 模式 `rational_precipitation`，用于拟合过饱和 gypsum 沉淀的 PHREEQC 累积体积曲线。
- 支持 `generation_realistic3` 和 `generation_realistic4` 的构造函数版 DFN 生成路径。
- 运行时会额外输出 `Output/DFN_raw.txt`，用于保存 backbone 之前的原始 DFN。

## 当前适用情况

当前版本最适合用于以下类型的问题：

- 过饱和 Ca-SO4 水在 gypsum fracture 中的 kinetic precipitation。
- 反应可以先用 PHREEQC 在一个参考水体积 `V_ref` 和参考反应面积 `A_ref` 上计算出累计矿物体积变化 `F_ref(t)`。
- 粒子追踪中的每个粒子只代表一部分注入水体积，使用 `V_particle / V_ref` 作为权重。
- 需要研究沉淀导致的 aperture reduction、流速改变、停留时间和空间分布之间的耦合。

当前版本不再建议用旧的“每个 particle parcel 单独生成一条 PHREEQC 曲线 `F_Vp(t)`”来解释反应。推荐做法是先生成一条固定的参考反应曲线 `F_ref(t)`，然后在 PT 中按粒子代表体积缩放。

## 技术栈

- 语言：C++
- 几何：CGAL
- 数值与特殊函数：Boost
- 线性代数/求解：uBLAS 风格结构与 SuiteSparse 相关依赖
- 可视化：OpenGL / GLUT
- 工程方式：Eclipse CDT + GNU Make

## 主流程

### 输入读取

`Input/File_names.txt` 当前主用为 3 行格式：

1. 域文件
2. 模拟文件
3. DFN 文件

程序仍兼容第 4 行 chemistry 文件名，但当前主路径不会再读取 chemistry input 文件。

`Parameters` 会进一步读取：

- `Input/Domain_files/` 中的域尺寸、基质扩散系数、孔隙度和左右边界水头；
- `Input/Simulation_files/` 中的粒子数、时间范围、反应时间步、输出时间步等；
- `Input/DFN_files/` 中的 DFN 模式与参数；
- `Input/Simulation_files/` 中末尾的 `Vref` 与 fracture thickness。

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

### 当前 chemistry-driven aperture 更新

当前主路径下，几何更新采用以下链条：

1. 每个粒子注入时被赋予 `representative_volume`；
2. 粒子在输运过程中累计 `particle_age`；
3. chemistry 模块计算参考系统的 `DeltaV_ref = F_ref(t_end) - F_ref(t_start)`；
4. 该体积增量按 `V_particle / V_ref` 缩放，当前 benchmark 关闭额外几何修正项；
5. 矿物体积变化换算为 segment 平均 aperture 增量；
6. aperture 更新后，若有裂隙关闭，则重建网络并重算流场。

当前 gypsum precipitation benchmark 使用的核心公式为：

```text
dV_particle = [F_ref(t_end) - F_ref(t_start)] * (V_particle / V_ref)
```

其中 `F_ref(t)` 来自 PHREEQC 的过饱和 gypsum 沉淀结果，并用 rational curve 表示：

```text
F_ref(t) = Pinf * [1 - w*(1 + a10*t)^(-1/9) - (1-w)*(1 + a2*t)^(-1)]
```

`Pinf` 为带符号矿物体积变化；在沉淀问题中 `Pinf < 0`，因此 aperture 会减小。aperture 的体积-几何转换公式没有为本 benchmark 单独改变。

旧的 `ChangeAperture()`、`ComputeAperture()` 等经验接口仍保留在代码里，但当前主路径默认不使用。

## V5-final benchmark

当前提交包含一个完整的过饱和 gypsum precipitation benchmark：

- 脚本：`Benchmark/gypsum_single_fracture_10cm_ap5mm/run_total_injection_reference_volume_precipitation_benchmark.py`
- 结果：`Output/benchmark_gypsum_single_fracture_10cm_ap5mm_total_injection_reference_volume_precipitation_20260426_174012/`

该 benchmark 使用 10 cm × 10 cm 单裂隙、初始 aperture 5 mm、固定两壁反应面积 `A_ref = 200 cm^2`。对不同 crossing time 计算总注入水体积 `V_ref = Q_in * t_injection`，每个 crossing-time case 只生成一条共享的 PHREEQC 参考曲线 `F_ref(t)`。输出目录中包含 PHREEQC 输入/输出、拟合曲线、PT 原始结果、aperture profile、aperture evolution、总沉淀量和 crossing-time 机理解释图。

## 核心文件

- `Code/src/PERFORM.cpp`
  - 主入口，负责参数读取、DFN 构建分流、backbone、输运和结果输出。
- `Code/src/Input_Output/Parameters.cpp`
  - 读取域、模拟和 DFN 参数，以及 simulation 中的 `Vref/thickness`。
- `Code/src/Domain_Definition/NetworkMeshes.cpp`
  - DFN 读取、DFN 生成、交点处理、backbone、流速计算。
- `Code/src/Transport/Transport.cpp`
  - 主输运循环、粒子注入、步进输运、chemistry 耦合 aperture 更新。
- `Code/src/Chemistry/Chemistry.cpp`
  - `DeltaV` 累积律、体积缩放和 aperture 增量换算。
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
