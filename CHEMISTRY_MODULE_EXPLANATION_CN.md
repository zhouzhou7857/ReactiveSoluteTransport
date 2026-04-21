# Chemistry 模块说明

## 1. 文件位置与职责

当前 chemistry 主路径的核心实现位于：

- [Chemistry.h](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/Chemistry/Chemistry.h)
- [Chemistry.cpp](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/Chemistry/Chemistry.cpp)

与其直接耦合的输运实现位于：

- [Particle.h](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/Transport/Particle.h)
- [Transport.cpp](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/Transport/Transport.cpp)

运行时参数读取与 chemistry 配置位于：

- [Parameters.h](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/Input_Output/Parameters.h)
- [Parameters.cpp](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/Input_Output/Parameters.cpp)
- [PERFORM.cpp](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/PERFORM.cpp)

## 2. 当前模块的设计目标

当前 chemistry 模块不再把主路径建立在“粒子浓度一阶衰减再换算 aperture”上，而是采用一个更直接的 PHREEQC 标定体积律：

1. 给定参考累计矿物体积变化 `DeltaV_ref(t)`
2. 用粒子的年龄区间 `t_start -> t_end` 计算本时间段的体积增量
3. 用 `V_particle / Vref` 把参考 parcel 的体积变化缩放到当前粒子
4. 可选再乘以几何相关修正因子
5. 把矿物体积变化转换为 segment aperture 变化

因此，这一模块当前更准确地说是：

`粒子年龄驱动的 DeltaV 累积律 + 体积缩放 + 几何反馈`

## 3. 当前实现的核心思想

### 3.1 粒子代表体积

在 [Particle.h](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/Transport/Particle.h) 中，每个粒子包含：

- `representative_volume`
- `particle_age`

其中：

- `representative_volume` 表示该数值粒子所代表的真实流体体积
- `particle_age` 表示该粒子已经经历的累计反应时间

### 3.2 参考累计体积律

当前活跃函数在 [Chemistry.cpp](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/Chemistry/Chemistry.cpp) 中由

- `GetCumulativeMineralVolumeChange(double t_particle)`

定义：

`DeltaV_ref(t) = A1 * (1 - exp(-k1 * t)) + A2 * (1 - exp(-k2 * t)) + L * t`

这条曲线表示：对一个参考水体积 `Vref`，当反应时间累积到 `t` 时，总矿物体积变化是多少。

### 3.3 单个粒子在一个时间段上的体积增量

当前主路径不是直接使用 `DeltaV_ref(t)`，而是使用其在相邻粒子年龄之间的差值：

`DeltaV_particle[t_start,t_end] = (DeltaV_ref(t_end) - DeltaV_ref(t_start)) * scale`

对应实现为：

- `ComputeParticleSegmentMineralVolumeChange(...)`

这里：

- `t_start = particle_age`
- `t_end = particle_age + dt_local`

也就是说，粒子每前进一段时间，就从累计曲线上取一段增量，而不是从头重复计算整条曲线。

## 4. 体积缩放如何实现

### 4.1 基本缩放

最基础的缩放是：

`scale = V_particle / Vref`

对应实现：

- `ComputeParticleVolumeScalingFactor(...)`

这表示：

- `PHREEQC` 或拟合公式给出的是参考水体积 `Vref` 的反应结果
- 当前粒子的实际贡献按它所代表的水体积线性缩放

### 4.2 可选修正项

当前代码还支持两个可选修正项：

1. effective diffusion-height factor
2. Vp-width correction

它们都通过环境变量在 `PERFORM.cpp` 中启用，不是默认主路径的一部分。

## 5. 从矿物体积变化到 aperture 变化

当前 aperture 更新不是通过旧的经验 `Nb/Nt` 公式完成，而是通过：

`delta_b = delta_V / (segment_length * thickness * 2)`

其中：

- `segment_length` 是当前 fracture segment 的长度
- `thickness` 是出平面厚度
- `/2` 对应双侧壁面均匀换算

实际实现位于：

- `ComputeApertureChangeFromMineralVolume(...)`

## 6. 粒子代表体积为什么是动态的

粒子不是固定携带常数体积，而是在注入时根据入口流动状态确定：

`V_particle = |u_inlet| * aperture_inlet * thickness * dt_particle`

实现位于：

- `ComputeParticleRepresentativeVolumeAtInjection(...)`

它表示：每个粒子代表该注入时间份额内通过入口裂缝的真实流体体积。

这使得在压力差边界条件下，粒子体积会随流速和 aperture 改变而动态调整。

## 7. 主路径中的几何更新逻辑

当前活跃逻辑位于 `Transport.cpp`：

1. `AccumulateParticleMineralVolumeChange(...)`
   - 根据粒子年龄区间和 `representative_volume` 计算 segment 上的 `DeltaV`
2. 每个 operator-splitting 时间步内，把所有粒子的 `DeltaV` 累计到对应 segment
3. `ApplyAccumulatedGeometryChange(...)`
   - 把 segment 总体积变化统一转换为 aperture 变化
4. 若某些 aperture 变为零，则触发网络重建与流场更新

## 8. 当前 chemistry 参数从哪里来

### 默认值

`A1 / k1 / A2 / k2 / L / Vref / thickness` 都在 `Chemistry.h` 中有默认值。

### 当前主用输入

当前主路径中：

- chemistry input 文件不再参与活跃路径
- `Vref` 和 `fracture thickness` 从 simulation input 读取
- 其他 `DeltaV` 系数默认来自代码

### 运行时覆盖

`PERFORM.cpp` 支持通过环境变量覆盖：

- `RST_CHEM_MODE`
- `RST_DELTA_V_A1`
- `RST_DELTA_V_K1`
- `RST_DELTA_V_A2`
- `RST_DELTA_V_K2`
- `RST_DELTA_V_L`
- `RST_DELTA_V_VREF`
- `RST_FRACTURE_THICKNESS`
- `RST_USE_EFFECTIVE_DIFFUSION_HEIGHT_FACTOR`
- `RST_EFFECTIVE_DIFFUSION_COEFFICIENT`
- `RST_EFFECTIVE_DIFFUSION_TIME`
- `RST_USE_VP_WIDTH_CORRECTION`

这套接口主要用于 benchmark 和敏感性分析。

## 9. 当前物理含义与适用解释

当前 chemistry 模块更适合解释为：

“每个粒子代表一份真实流体体积；该流体包在裂缝中随着粒子年龄增长，按参考 `DeltaV` 累积律产生矿物体积变化；这些体积变化再反馈到裂缝 aperture。”

它更接近：

- PHREEQC 标定的等效反应 parcel 模型
- 粒子年龄驱动的累计反应模型
- 体积守恒驱动的几何演化模型

而不是：

- 显式多组分平衡-动力学求解器
- 当前主路径上的浓度衰减控制模型

## 10. 已实现部分

明确已实现：

1. 参考累计体积律 `DeltaV_ref(t)`
2. 粒子年龄区间增量计算
3. `V_particle / Vref` 体积缩放
4. 可选几何修正因子
5. segment 上累计矿物体积变化
6. 体积变化到 aperture 变化的换算
7. aperture 更新后反馈到流场与输运

## 11. 保留但非主路径的旧功能

以下内容仍存在，但不是当前默认主路径：

1. `UpdateReactiveConcentration()`
2. `EvaluateReactiveStep()`
3. 基于 `Nb/Nt` 的经验 aperture 更新
4. chemistry input 文件驱动的旧工作流

因此，阅读当前代码时，应优先以 `ComputeParticleSegmentMineralVolumeChange()` 和 `AccumulateParticleMineralVolumeChange()` 为主线理解 chemistry 模块。
