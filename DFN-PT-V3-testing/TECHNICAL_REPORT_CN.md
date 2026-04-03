# 基于实现的技术报告

## 说明

本报告严格基于当前有效代码 `Code/src`。除非用于说明哪些实现已不在主路径上，否则不讨论备份目录、历史版本和结果后处理脚本。

## A. 代码整体目标

该代码实现了二维离散裂隙网络（DFN）中的粒子追踪溶质运移，并将输运过程与 DFN 结构演化耦合。当前主线实现中，结构演化不是由通用化学平衡或浓度场直接驱动，而是由“粒子穿越裂隙网格的有效通量指标”驱动裂隙开度更新，再基于新开度重算流场。

当前代码中最核心的耦合链条是：

- 粒子输运生成每个裂隙网格的粒子穿越量；
- 粒子穿越量驱动局部裂隙开度更新；
- 开度变化改变水力传导能力和等效传输率；
- DFN 全局流场重新求解；
- 下一时间步的粒子输运使用更新后的 DFN。

## B. 程序主流程

### 主入口文件、核心模块与关键函数

主入口与总控：

- `Code/src/PERFORM.cpp`

核心模块：

- `Code/src/Input_Output/Parameters.h/.cpp`
- `Code/src/Domain_Definition/Domain.h/.cpp`
- `Code/src/Domain_Definition/NetworkMeshes.h/.cpp`
- `Code/src/Domain_Definition/DFNComputation.h/.cpp`
- `Code/src/Domain_Definition/FractureMesh.h/.cpp`
- `Code/src/Transport/Transport.h/.cpp`
- `Code/src/Transport/Projection.cpp`
- `Code/src/Transport/Transfer.cpp`
- `Code/src/Chemistry/Chemistry.h/.cpp`
- `Code/src/Input_Output/Results.h/.cpp`

辅助模块：

- `Code/src/Visualisation/*`
- `Code/src/RngStream/*`
- `Code/src/Utilitaries/*`

### 从输入到输出的执行流程

1. `PERFORM.cpp` 中的 `main()` 创建 `Parameters param` 和 `Domain domain`。
2. `Parameters::Parameters()` 先读取 `Input/File_names.txt`，再调用 `read_param()`。
3. `Parameters::read_param()` 读取：
   - 计算域尺寸、基质扩散系数、孔隙度、左右边界水头；
   - 粒子数、转移概率控制、矩阵模式、突破曲线时间范围、随机种子、注入持续时间、输出间隔、反应时间步；
   - DFN 文件中的生成/读取模式。
4. `NetworkMeshes` 根据选定 DFN 文件构建裂隙网络。
5. 如果 `backbone=true`，则调用 `return_backbone()`：
   - 剔除不连通部分；
   - 基于流速阈值保留骨架；
   - 内部调用 `ComputeFlowVelocities()` 求初始流场。
6. 构造 `Transport transport(...)`。
7. 调用 `transport.Particles_Transport(arrival_times, option_injection)` 执行主输运-结构演化循环。
8. 若输运成功：
   - `Results::post_processing()` 生成突破曲线统计；
   - `Results::writing()` 输出 `cdf.txt` 和 `pdf.txt`；
   - `WritePositionSnapshotsCSV()` 输出粒子位置快照。
9. 最后写出：
   - 当前 DFN；
   - 初始 DFN；
   - 开度变化后的 DFN 比较文件。

### 活跃路径与非主路径

当前主程序实际使用的是：

- `Particles_Transport()`
- `particle_displacement_step()`

旧函数如 `particle_displacement()` 仍保留在代码中，也包含开度更新逻辑，但它不属于当前 `PERFORM.cpp` 的主执行路径。

## C. 粒子追踪算法：实现步骤

### 1. 粒子初始化

粒子状态由 `Particle` 类定义，主要变量包括：

- 当前位置 `M`
- 当前裂隙网格编号 `mesh_index`
- 前一网格编号 `prev_mesh_index`
- 网格访问历史
- 当前时间 `t`
- 注入时间 `t_injection`
- 当前裂隙内累计停留时间 `t_in_fract`
- 当前裂隙内累计平流距离 `L_in_fract`

初始化发生在 `Transport::Particles_Injection()` 中。

已实现行为：

- 创建 `nb_part` 个粒子；
- 初始空间位置并不立即设置，先置为未定义；
- 若 `t_injection > 0`，粒子注入时间在 `[0, t_injection]` 上均匀展开；
- 若 `t_injection <= 0`，所有粒子时间初始为 `0`。

因此当前实现是“先初始化时间，再在主循环中按时刻注入空间位置”。

### 2. 粒子注入

粒子注入在 `Particles_Transport()` 中真正完成。

首先由 `CollectInputPositions()` 从当前 DFN 中识别入口裂隙。入口条件是裂隙某个端点位于左边界。

当前实现支持三种注入方式：

- `option_injection = 0`：左边界全入口均匀注入；
- `option_injection = 1`：选择距离左边界中点最近的单一入口进行局部注入；
- `option_injection = 2`：按入口裂隙的 `aperture * |velocity|` 进行通量加权注入。

`ComputeInjectionCounts()` 负责把总粒子数分配到各入口裂隙：

- 先按权重计算期望值；
- 取整数部分；
- 剩余粒子按最大余数分配。

当粒子满足：

- `mesh_index == -1`
- 且 `t <= t_current`

时，该粒子在当前反应步被真正投放到某个入口位置和入口裂隙。

### 3. 粒子在裂隙中的运动

当前活跃的单步推进函数是 `particle_displacement_step()`。它将单个粒子推进到目标时间：

- `t_target = pa.t + dt_step`

代码实现了两种输运模式：

- 无限基质模式：`infinite_matrix_displacement_step()`
- 有限基质模式：`finite_matrix_displacement_step()`

#### 无限基质模式

已实现步骤：

1. 用 `distance_from_M_to_extremity()` 找到顺流方向终点。
2. 用 `Mesh_Time_Advection()` 计算平流时间。
3. 用 `Get_Total_Time_From_Advec_Time1D()` 将平流时间转换为包含基质扩散效应的总停留时间。
4. 若该事件会超过 `t_target`，则只在裂隙内推进到对应的截断位置。
5. 否则粒子到达当前裂隙末端。

当前实现的总停留时间公式写在 `FractureMesh::Get_Total_Time_From_Advec_Time1D()` 中：

- `B = sqrt(Dm) * porosity * advection_time / (0.5 * aperture)`
- `res_time = advection_time + 0.25 * (B / erfc_inv(u))^2`

其中 `u` 为随机数。

这说明当前主路径中的停留时间修正明确使用了 1D 基质扩散形式。

#### 有限基质模式

已实现步骤：

1. 确定顺流方向的裂隙终点 `M_out`。
2. 用 `Orthogonal_Projection_M()` 构造从当前裂隙到邻近裂隙的正交投影。
3. 用 `Advection_Time_Computation()` 根据局部转移概率上限 `proba_transfer` 求稳定平流子步长。
4. 用 `Get_Total_Time_From_Advec_Time1D()` 得到总停留时间。
5. 将超出平流部分解释为扩散时间 `t_diff`。
6. 用 `Transfer_Probability_and_Transfer_Time()` 判断在 `t_diff` 内是否发生跨裂隙转移。
7. 若发生转移：
   - 确定转移时间；
   - 确定到达哪一个屏障/裂隙；
   - 将粒子投影到目标裂隙位置。
8. 若不发生转移，则粒子留在原裂隙继续运动。
9. 若整个事件会超过 `t_target`，则在当前时间步截断。

因此，有限基质模式不仅在节点处切换路径，还实现了通过基质向邻近裂隙扩散转移的机制。

### 4. 路径选择

当粒子到达裂隙端点且尚未离开域时，下一条裂隙由 `Mesh_Neighbour_Choice_And_Test()` 选择。

实现过程如下：

1. 判断粒子当前更接近当前裂隙的哪个端点。
2. 收集与该节点相连的候选邻接裂隙。
3. 对每条候选裂隙计算“流量代理”：
   - 若流向从该节点流出，则用 `velocity * aperture`；
   - 若该裂隙方向相反，则用 `-velocity * aperture`。
4. 调用 `select_neighbour_slr()` 将这些值转换成概率。

`select_neighbour_slr()` 实现的是一种 streamline routing 规则：

- 若无有效下游路径：
  - 若已在右边界，则可视为出口；
  - 否则报警；
- 若只有一条路径：概率为 1；
- 若有两条路径：
  - 若存在“直行”路径，则优先比较直行与分叉；
  - 否则按流量比分配；
- 若三条或以上路径：按流量比分配。

因此当前实现并不是最短路径或纯几何路径选择，而是基于流量的路径选择。

### 5. 交点处理

当前代码中有两类“交点/切换”处理：

- 基于 DFN 拓扑节点的裂隙切换：`Mesh_Neighbour_Choice_And_Test()`
- 基于基质扩散的近邻裂隙转移：`Orthogonal_Projection_M()` + `Transfer_*`

所以实现上同时支持：

- 沿裂隙网络图结构传播；
- 通过基质扩散跨越到邻近但不一定直接相连的裂隙。

### 6. 行程时间与停留时间计算

当前实现中与时间相关的核心计算包括：

- 裂隙内平流时间：`Mesh_Time_Advection()`
- 包含扩散修正的总停留时间：`Get_Total_Time_From_Advec_Time1D()`
- 有限基质转移的概率与时间：
  - `Transfer_Probability_Computation()`
  - `Transfer_Time_Computation()`
  - `Transfer_Time_FPTD()`
  - `Transfer_Time_Feller()`

从 `Transfer.cpp` 可见：

- 若只有一个屏障，使用 first-passage time 分布；
- 若有两个屏障，使用 Feller 两屏障分布。

### 7. 边界处理

当前已明确实现的边界逻辑：

- 左边界：入口识别与粒子注入；
- 右边界：粒子离开系统的判据；
- 上下边界：流动求解中默认采用零 Neumann 条件。

`Projection.cpp` 中存在 `ReflectionOnBorder()` 等反射辅助函数，说明代码结构上考虑过顶/底边界反射扩散。但在当前主执行路径 `finite_matrix_displacement_step()` 中，并不能明确看到这些反射函数被直接调用。

因此可以明确区分：

- 顶/底边界反射能力在辅助函数层面存在；
- 但它是否在当前活跃主路径中真正起作用，并不完全明确。

### 8. 粒子退出与输出记录

粒子退出逻辑位于 `Particles_Transport()`：

- 若粒子到达右边界，则记录 `arrival_times[particle_id] = arrival_time`；
- 然后将粒子标记为失活。

输出记录包括：

- `RecordPositions()`：把粒子快照存入内存；
- `WriteDFNSnapshot()`：写出 `DFN_step*.txt`；
- `WritePositionSnapshotsCSV()`：写出 `particle_positions_t*.csv`；
- `Results::post_processing()`：生成突破曲线统计；
- `Results::writing()`：写出 `cdf.txt` 和 `pdf.txt`。

## D. 关键粒子相关变量及其物理含义

- `nb_part`
  - 总注入粒子数。
- `proba_transfer`
  - 有限基质模式中用于控制转移概率上限的参数。
- `simu_option`
  - 输运模式选择：
    - `0`：无限基质
    - `1`：有限基质
- `t_injection`
  - 粒子注入持续时间。
- `output_interval`
  - 粒子位置和 DFN 快照输出间隔。
- `reaction_dt`
  - DFN 更新与流场重算的固定耦合时间步。
- `pa.t`
  - 粒子当前时间。
- `pa.t_injection`
  - 粒子的注入时刻。
- `pa.t_in_fract`
  - 粒子自进入当前裂隙以来在该裂隙中的累计停留时间。
- `pa.L_in_fract`
  - 粒子自进入当前裂隙以来在该裂隙中的累计平流距离。
- `moved_distances[mesh_id]`
  - 某一反应步内粒子在某裂隙中的累计平流距离。
- `step_full_count[mesh_id]`
  - 当前反应步内对某裂隙完成的整段穿越次数。
- `step_partial_sum[mesh_id]`
  - 当前反应步内对某裂隙的分数穿越量累计。
- `Nb_eff`
  - 某裂隙在一个反应步内的有效粒子穿越量。
- `reference_injected_count_dt`
  - 一个反应时间步内参考注入粒子数，用于归一化开度更新。

## E. 粒子输运如何影响 DFN 结构

### 明确实现的内容

#### 裂隙开度

已实现：是。

更新位置：

- `Transport::Particles_Transport()` 计算每个裂隙的 `Nb_eff`；
- `NetworkMeshes::ChangeAperture()` 执行开度更新；
- `Chemistry::ComputeAperture()` 计算新开度。

更新公式：

- `ratio = min(1, Nb / Nt)`
- `delta = K_REACTION * dt * ratio`
- `b_new = b_old - delta`

其中 `K_REACTION = -1.5e-10`。

由于该常数为负，当前公式在 `ratio > 0` 时会使开度增大。这与代码注释中“溶解导致开度增大”的意图一致。

更新特征：

- 局部：逐裂隙网格更新；
- 周期性：每个 `reaction_dt` 更新一次；
- 累积性：在整个模拟中持续累计；
- 阈值敏感：开度被限制为非负，后续会检测是否为零。

#### 孔隙度

已实现动态更新：否。

证据：

- `porosity` 从输入文件读入；
- 被输运时间与扩散公式使用；
- 之后没有任何代码修改它。

#### 渗透率

作为显式动态状态量更新：否。

代码中没有单独的可演化渗透率字段。

#### 传输率（transmissivity）

作为独立持久状态显式更新：否。

作为开度函数隐式变化：是。

证据：

- `HydraulicProperties::ReturnTransmissivity()` 按立方律使用 `aperture^3`；
- 每次流场重算时，新的开度会自动影响 transmissivity；
- 但代码没有独立维护一个随时间演化的 transmissivity 场。

#### 局部几何

坐标或裂隙线段几何形状更新：否。

裂隙端点坐标不随反应更新。变化的是开度，不是线段本身的几何位置。

#### 裂隙段属性

部分已实现：是。

会变化的属性：

- aperture；
- velocity（在流场重算后）；
- 有效连通状态（在剔除零开度裂隙后）。

不会变化的属性：

- 端点坐标；
- 裂隙长度；
- 基质孔隙度。

#### DFN 连通性与拓扑

已实现：是，但属于事件触发型。

机制：

- 若某裂隙开度为零，则视为关闭；
- 零开度裂隙会被移除；
- 网络重建；
- 节点重新编号；
- 重新提取 backbone；
- 已存在粒子可能被重映射或直接移除。

主要位置：

- `RebuildNetworkRemovingZeroAperture()`
- `return_backbone()`
- `RebuildNodeNumbering()`
- `Particles_Transport()` 中的粒子重映射逻辑

这不是连续几何演化，而是阈值触发的拓扑裁剪。

### 仅由代码结构暗示但未完全展开的内容

- 代码注释暗示其物理意图是反应溶解导致裂隙增宽，再改变流动与输运。
- 但当前实现本质上是经验型、粒子通量驱动的开度模型，而不是由质量守恒、浓度场或反应动力学方程直接推导。

### 明确缺失的内容

- 没有孔隙度演化；
- 没有矿物质量守恒；
- 没有浓度场求解；
- 没有显式渗透率张量；
- 没有裂隙几何坐标的演化；
- 没有力学变形；
- 没有完整的沉淀/堵塞模型。

## F. 输运与 DFN 演化之间的反馈

当前代码确实实现了双向反馈。

### 正向：输运影响 DFN

在每个反应步结束时：

1. 统计每个粒子在各裂隙中的平流距离；
2. 将其转换为每个裂隙的有效粒子穿越量 `Nb_eff`；
3. 用 `Nb_eff / reference_injected_count_dt` 更新局部开度。

### 反向：DFN 影响输运

在开度更新后：

1. `ComputeFlowVelocities()` 在更新后的 DFN 上重算流速；
2. `NormalizeMeshDirections()` 将流向统一为正向编号；
3. 若裂隙关闭或流速极低，则触发网络重建；
4. `net_mesh = net_mesh_modified`，使下一反应步使用更新后的 DFN；
5. 因此下一步的粒子输运、入口权重和路径选择都依赖新流场。

### 耦合性质

- 开度更新：局部、周期性；
- 流场重算：全局、周期性；
- 拓扑变化：局部触发、全局重建；
- 时间推进方式：显式分裂耦合，而非连续同步耦合。

换言之，粒子在一个 `reaction_dt` 内是在“冻结 DFN”上运动；DFN 只在步间更新。

## G. 主要假设与简化

- 裂隙被表示为二维线段，带标量开度。
- 流动仅在 DFN 图结构上求解。
- 基质效应通过扩散修正停留时间和跨裂隙转移体现，而不是通过显式基质浓度场体现。
- 开度演化由经验公式驱动，不由溶质浓度或反应方程直接决定。
- 基质孔隙度和扩散系数为空间常数。
- DFN 只在固定 `reaction_dt` 上更新。
- 节点处分流采用基于流量的 streamline routing 规则。

## H. 当前局限与缺失部分

- 耦合变量是粒子穿越量，而不是溶质质量或浓度。
- 没有显式化学组分、反应物或平衡求解。
- 没有孔隙度更新。
- 没有独立的渗透率演化变量。
- 没有连续的几何演化，拓扑变化只在阈值触发时发生。
- 没有裂隙坐标增宽或形态变化。
- 当前有限基质主路径仍调用 `Get_Total_Time_From_Advec_Time1D()`，说明停留时间修正仍沿用 1D 扩散形式。
- 代码中存在一些边界处理辅助函数，但它们是否都在当前活跃路径中起作用并不完全明确。
- 源码中保留了旧版输运函数和注释掉的替代方案，说明实现历史仍混在当前源码中。
- 活跃代码存在潜在健壮性问题，例如 `fabs(current_mesh.velocity==0.0)` 看起来是对布尔表达式取绝对值，而不是对速度本身取绝对值。

## I. 文本流程图

```text
Input/File_names.txt
    ->
Parameters::read_param()
    ->
读取计算域、模拟参数、DFN 文件
    ->
构建 NetworkMeshes
    ->
return_backbone()
    ->
ComputeFlowVelocities() 得到初始流场
    ->
Transport::Particles_Transport()
    ->
Particles_Injection()
    ->
对每个反应时间步:
    ->
将达到注入时刻的粒子投放到左边界入口裂隙
    ->
对每个活跃粒子:
    ->
particle_displacement_step()
    ->
    infinite_matrix_displacement_step()
    或
    finite_matrix_displacement_step()
    ->
    若到达节点:
        Mesh_Neighbour_Choice_And_Test()
    ->
    若到达右边界:
        记录到达时间并失活
    ->
累积各裂隙的 moved_distances
    ->
计算每个裂隙的 Nb_eff
    ->
ChangeAperture(mesh, reaction_dt, Nb_eff, reference_injected_count_dt)
    ->
若出现零开度或极低流速:
    重建 DFN 与 backbone
    重映射或移除粒子
    ->
ComputeFlowVelocities() 重算更新后 DFN 的流场
    ->
进入下一反应步
    ->
Results::post_processing()
    ->
Results::writing()
    ->
WritePositionSnapshotsCSV() 与 DFN 输出
```

## J. 简短要点总结

- 当前主线是一个带周期性开度更新的 DFN 粒子追踪模型。
- 粒子在每个固定反应步内在冻结 DFN 上运动。
- 每步结束后，粒子穿越量被转换为每个裂隙的有效穿越指标 `Nb_eff`。
- `Nb_eff` 局部更新开度。
- 开度变化通过流场重算反馈到后续粒子输运。
- 若裂隙关闭或失去水力活性，DFN 会被重建，粒子可能被重映射或移除。
- 当前没有动态孔隙度、浓度场或完整反应化学模型。

## 理解该模型最重要的源文件

- `Code/src/PERFORM.cpp`
- `Code/src/Transport/Transport.cpp`
- `Code/src/Transport/Transfer.cpp`
- `Code/src/Transport/Projection.cpp`
- `Code/src/Domain_Definition/NetworkMeshes.cpp`
- `Code/src/Domain_Definition/DFNComputation.cpp`
- `Code/src/Domain_Definition/FractureMesh.cpp`
- `Code/src/Chemistry/Chemistry.cpp`
- `Code/src/Input_Output/Parameters.cpp`
- `Code/src/Input_Output/Results.cpp`

## 粒子输运的主要函数/子程序

- `Transport::Particles_Transport()`
- `Transport::Particles_Injection()`
- `Transport::particle_displacement_step()`
- `Transport::infinite_matrix_displacement_step()`
- `Transport::finite_matrix_displacement_step()`
- `Transport::Mesh_Neighbour_Choice_And_Test()`
- `Transport::select_neighbour_slr()`
- `Transport::Transfer_Probability_and_Transfer_Time()`
- `Transport::Transfer_Probability_Computation()`
- `Transport::Transfer_Time_Computation()`
- `FractureMesh::Get_Total_Time_From_Advec_Time1D()`
- `FractureMesh::Advection_Time_Computation()`

## DFN 更新的主要函数/子程序

- `NetworkMeshes::ChangeAperture()`
- `Chemistry::ComputeAperture()`
- `RebuildNetworkRemovingZeroAperture()`
- `NetworkMeshes::return_backbone()`
- `RebuildNodeNumbering()`
- `ComputeFlowVelocities()`
- `NetworkMeshes::EvaluateFlowVelocities()`
- `ZeroApertureIfLowVelocity()`

## 当前耦合策略的简短总结

当前耦合策略属于算子分裂、粒子穿越量驱动的 DFN 演化方法。在每个固定反应时间步内，粒子在冻结流场和冻结 DFN 上运动；该时间步结束后，代码根据粒子在各裂隙中的有效穿越量局部更新开度，再用更新后的开度全局重算 DFN 流场。若某些裂隙关闭或失去水力活性，则进一步重建网络并在下一输运步中使用新 DFN。
