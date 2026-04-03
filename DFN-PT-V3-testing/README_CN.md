# ReactiveSoluteTransport

## 项目概述

ReactiveSoluteTransport 是一个基于离散裂隙网络（DFN）的二维溶质运移数值模拟项目，主程序位于 `Code/src`。程序从输入文件读取计算域、裂隙网络和模拟参数，先求解裂隙网络中的流动，再执行粒子追踪，并在时间步推进过程中按照粒子穿越强度更新裂隙开度，必要时重建网络连通性与流场。

当前代码的核心特点如下：

- 以粒子追踪方式模拟裂隙中的平流和基质扩散。
- 支持无限基质与有限基质两种扩散选项。
- 支持连续注入、局部注入和入口通量加权注入。
- 在输运过程中按固定反应时间步更新裂隙开度，并在开度归零后重建 DFN 骨架。
- 输出突破曲线结果，同时记录粒子位置快照和 DFN 开度演化快照，便于后处理对比。

## 技术栈

### 前端

- 无独立前端应用。
- `Code/src/Visualisation` 提供基于 OpenGL/GLUT 的本地可视化代码，用于裂隙网络和结果展示。

### 后端

- C++ 数值模拟程序。
- 基于类与头文件接口组织模块，没有 Web 服务或 RPC 接口。
- 使用 Boost 的特殊函数支持部分输运概率和时间分布计算。
- 使用 CGAL 几何类型与几何运算处理点、线段、投影和边界关系。
- 使用 uBLAS 风格矩阵/向量封装与线性系统求解工具处理流动计算。

### 数据与存储

- 无数据库。
- 输入数据以纯文本文件组织在 `Input/` 下。
- 输出结果以 `txt`、`csv`、`png`、`tsv`、`xlsx` 等文件形式写入 `Output/`。

### 开发工具

- Eclipse CDT 工程配置文件：`Code/.project`、`Code/.cproject`。
- GNU 工具链配置痕迹明显，工程配置中包含 GCC/G++、GNU Make。
- 仓库中保留 `Debug/`、`Release/`、`.vscode/` 等本地开发环境产物。

## 项目架构

### 整体架构

项目采用单体式数值模拟架构，入口程序负责顺序编排“参数读取 -> 几何/DFN 构建 -> 流动计算 -> 粒子输运 -> 化学开度更新 -> 结果输出”。

核心模块关系如下：

- `Input_Output` 读取域、模拟与边界参数，并负责突破曲线结果写出。
- `Domain_Definition` 定义计算域、裂隙网格、节点编号、邻接关系和 DFN 流动求解所需结构。
- `Flow` 与 `DFNComputation` 负责裂隙网络水头和流速计算。
- `Transport` 负责粒子注入、步进输运、跨裂隙转移、粒子轨迹记录和反应步调度。
- `Chemistry` 根据反应规则计算裂隙开度变化。
- `Utilitaries` 提供几何、随机数、线性代数、常量和通用数值工具。
- `Visualisation` 提供可选的网络与结果显示代码。

### 数据流向

1. `Code/src/PERFORM.cpp` 创建 `Parameters` 和 `Domain`，通过 `Input/File_names.txt` 选择当前案例文件。
2. `Parameters` 从 `Input/Domain_files`、`Input/Simulation_files`、`Input/DFN_files` 读取域大小、扩散参数、时间步、注入时长和 DFN 来源。
3. `NetworkMeshes` 读取或生成裂隙网络，构建节点、边界和裂隙映射。
4. 程序通过 `return_backbone` 和 `ComputeFlowVelocities` 提取连通骨架并计算正向流速。
5. `Transport::Particles_Transport` 以反应时间步推进粒子，统计各裂隙穿越量，并调用 `NetworkMeshes::ChangeAperture` 更新开度。
6. 当某些裂隙开度降为零时，代码会移除失效网格、重建节点编号、重新提取骨架并重新计算流场。
7. `Results` 输出 `cdf.txt` 与 `pdf.txt`，`Transport` 输出粒子位置快照和 DFN 时间快照。

## 目录结构

以下目录树聚焦当前主线代码和主要数据目录，省略 `.git`、`.venv`、`node_modules`、`dist` 等无关目录：

```text
ReactiveSoluteTransport/
├── README.md
├── Code_Update_Log.txt
├── test_case_matrix.xlsx
├── Code/
│   ├── .project
│   ├── .cproject
│   ├── src/
│   │   ├── PERFORM.cpp
│   │   ├── Chemistry/
│   │   │   ├── Chemistry.h
│   │   │   └── Chemistry.cpp
│   │   ├── Domain_Definition/
│   │   │   ├── Domain.h/.cpp
│   │   │   ├── FractureMesh.h/.cpp
│   │   │   ├── HydraulicProperties.h/.cpp
│   │   │   ├── NetworkMeshes.h/.cpp
│   │   │   └── DFNComputation.h/.cpp
│   │   ├── Flow/
│   │   │   └── FlowComputation.h/.cpp
│   │   ├── Input_Output/
│   │   │   ├── Parameters.h/.cpp
│   │   │   ├── BoundaryConditions.h/.cpp
│   │   │   └── Results.h/.cpp
│   │   ├── Transport/
│   │   │   ├── Transport.h/.cpp
│   │   │   ├── Particle.h/.cpp
│   │   │   ├── Projection.cpp
│   │   │   └── Transfer.cpp
│   │   ├── Utilitaries/
│   │   │   ├── Constantes.h
│   │   │   ├── LinearSystem.h/.cpp
│   │   │   ├── RandomNumber.h/.cpp
│   │   │   ├── Segment.h/.cpp
│   │   │   ├── Point_Cgal.h/.cpp
│   │   │   ├── FluxPoint2D.h/.cpp
│   │   │   ├── Structures.h/.cpp
│   │   │   ├── UblasStructures.h/.cpp
│   │   │   └── 其他几何与数值工具
│   │   ├── Visualisation/
│   │   │   ├── DFNVisu.h/.cpp
│   │   │   └── DisplayResults.h/.cpp
│   │   └── RngStream/
│   │       └── rngstream.h/.cpp
│   ├── Backup0318/
│   │   └── src/
│   ├── Backup0318_unmodified/
│   │   └── src/
│   ├── Debug/
│   ├── Release/
│   └── input/
├── Input/
│   ├── File_names.txt
│   ├── Domain_files/
│   │   ├── Domain1.txt
│   │   ├── Domain2.txt
│   │   ├── Domain3.txt
│   │   ├── Domain_gypsum.txt
│   │   └── Domain_Description.txt
│   ├── Simulation_files/
│   │   ├── Simulation.txt
│   │   ├── Simulation_gypsum.txt
│   │   ├── Simulation_Description.txt
│   │   └── Test-0324/
│   └── DFN_files/
│       ├── DFN_newformat_17.txt
│       ├── DFN.txt
│       ├── DFN_gypsum*.txt
│       ├── fracture_data_*.txt
│       ├── gypsum discretization/
│       └── DFN_Description.txt
└── Output/
    ├── test/
    ├── Compare/
    ├── InitialResults/
    ├── RunSimulation 1/
    ├── 0324/
    └── 其他历史结果目录
```

目录用途说明：

- `Code/src`：当前主线实现，README 中的功能说明均以这里为准。
- `Code/Backup0318`、`Code/Backup0318_unmodified`：历史快照，不应与当前运行逻辑混用。
- `Input`：案例选择文件和各类输入参数文件。
- `Output`：模拟结果、比较图和实验性后处理脚本输出目录。

## 核心文件说明

### 项目入口文件和配置文件

- `Code/src/PERFORM.cpp`
  - 项目主入口。
  - 负责创建参数对象、读取 DFN、执行骨架提取、启动粒子输运，并在结束时写出 DFN 与统计结果。
- `Input/File_names.txt`
  - 当前案例选择器。
  - 按顺序指定域文件、模拟文件和 DFN 文件名，`Parameters` 构造函数会直接读取它。
- `Input/Domain_files/Domain1.txt`
  - 当前默认案例的计算域与基质参数文件。
  - 包含域尺寸、基质扩散系数、孔隙度和左右边界水头。
- `Input/Simulation_files/Simulation.txt`
  - 当前默认案例的输运与反应调度参数。
  - 包含粒子数、裂隙转移概率、矩阵扩散模式、输出时间范围、随机种子、注入持续时间、输出间隔和反应时间步。
- `Input/DFN_files/DFN_newformat_17.txt`
  - 当前默认案例的 DFN 输入。
  - 文件首行声明 `file`，第二行给出裂隙数量，后续各行给出裂隙两端坐标和开度。
- `Code/.project`、`Code/.cproject`
  - Eclipse CDT 工程配置。
  - 从配置内容可见工程依赖 CGAL、Boost、SuiteSparse/UMFPACK、OpenGL/GLUT 等本地库。

### 核心业务逻辑实现

- `Code/src/Transport/Transport.cpp`
  - 当前最核心的业务实现文件。
  - `Particles_Transport` 负责粒子注入、逐反应时间步推进、出口到达判定、裂隙穿越统计、开度更新、失效裂隙剔除、骨架重建和流场重算。
  - 同文件还负责位置快照与 DFN 时间快照输出。
- `Code/src/Transport/Projection.cpp`
  - 处理裂隙间正交投影和周期边界相关投影逻辑，是粒子跨裂隙转移的重要几何基础。
- `Code/src/Transport/Transfer.cpp`
  - 计算粒子通过基质向邻近裂隙转移的概率与转移时间分布。
- `Code/src/Domain_Definition/NetworkMeshes.cpp`
  - 负责从文件构建 DFN、裂隙插入、节点编号、交点计算、骨架提取，以及开度变化后的网格更新。
  - `ChangeAperture` 和 `ChangeApertureByRatio` 直接作用于裂隙开度。
- `Code/src/Domain_Definition/DFNComputation.cpp`
  - 负责根据边界条件组装并求解 DFN 流动问题，计算水头和流速。
- `Code/src/Domain_Definition/FractureMesh.cpp`
  - 封装单条裂隙网格的时间、距离、流速和扩散相关计算。

### 数据模型和 API 接口

项目没有对外的 HTTP API，核心“接口”主要是 C++ 类与函数接口：

- `Code/src/Input_Output/Parameters.h`
  - `Parameters` 是全局配置数据模型，统一承载域、输运、DFN 和边界条件参数。
- `Code/src/Transport/Particle.h`
  - `Particle` 是粒子状态模型，记录当前位置、当前裂隙、历史轨迹、注入时间和在裂隙中的累计停留时间。
- `Code/src/Domain_Definition/FractureMesh.h`
  - `FractureMesh` 是单条裂隙段的数据模型，包含端点、流速、开度、邻接网格和原始编号。
- `Code/src/Domain_Definition/NetworkMeshes.h`
  - `NetworkMeshes` 是 DFN 主数据模型，维护裂隙集合、节点映射、边界映射、裂隙索引和流动相关状态。
- `Code/src/Input_Output/BoundaryConditions.h`
  - `BoundaryConditions`、`BoundaryConditionsDFN`、`BoundaryConditionsDef` 定义边界条件的数据接口。
- `Code/src/Input_Output/Results.h`
  - `Results` 提供突破曲线后处理与写出接口。

### 关键组件和服务模块

- `Code/src/Chemistry/Chemistry.h`
  - 定义反应常数 `K_REACTION` 以及 `ComputeAperture`、`ComputeDeltaAperture`、`safe_ratio`。
  - 当前实现按每个反应时间步的有效粒子通量比例更新裂隙开度，并将结果钳制为非负。
- `Code/src/Input_Output/Results.cpp`
  - 将粒子到达时间转换为累计分布 `cdf.txt` 和概率密度 `pdf.txt`。
- `Code/src/Utilitaries/LinearSystem.cpp`
  - 为 DFN 流动求解提供线性系统计算能力。
- `Code/src/RngStream/rngstream.cpp`
  - 提供随机数流，供粒子扩散和随机转移过程使用。
- `Code/src/Visualisation/DFNVisu.cpp` 与 `Code/src/Visualisation/DisplayResults.cpp`
  - 提供本地图形展示能力，但在当前主流程中属于辅助模块，不是必须路径。

## 当前默认案例补充说明

根据当前仓库输入文件：

- 默认域文件是 `Domain1.txt`，域尺寸为 `0.00746 m x 0.00586 m`。
- 默认模拟文件是 `Simulation.txt`，当前配置为 `50000` 个粒子，反应时间步为 `100 s`，位置输出间隔为 `1e5 s`。
- 默认 DFN 文件是 `DFN_newformat_17.txt`，当前文件声明共有 `568` 条裂隙记录。

## 输出结果说明

当前主流程会在 `Output/` 下生成或覆盖以下典型结果：

- `cdf.txt`、`pdf.txt`：突破曲线统计结果。
- `DFN.txt`、`DFN_init.txt`、`DFN_aperture_delta.txt`：当前 DFN、初始 DFN 和开度变化结果。
- `DFN_step*.txt`：反应时间步对应的 DFN 演化快照。
- `particle_positions_t*.csv`：按输出时刻记录的粒子位置快照。

这些输出已经在仓库中的 `Output/test`、`Output/RunSimulation 1`、`Output/Compare` 等目录留下样例。
