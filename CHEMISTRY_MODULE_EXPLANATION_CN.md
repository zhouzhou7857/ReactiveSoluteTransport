# 化学反应模块说明

## 1. 文件位置与职责

当前代码中，化学反应相关的核心实现位于：

- [Chemistry.h](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/Chemistry/Chemistry.h)
- [Chemistry.cpp](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/Chemistry/Chemistry.cpp)

与其直接耦合的输运实现位于：

- [Particle.h](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/Transport/Particle.h)
- [Transport.cpp](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/Transport/Transport.cpp)

输入读取与运行时 chemistry 参数配置位于：

- [Parameters.h](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/Input_Output/Parameters.h)
- [Parameters.cpp](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/Input_Output/Parameters.cpp)
- [PERFORM.cpp](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/PERFORM.cpp)

## 2. 模块的设计目标

当前化学模块的目标不是构建一个完整的多组分地球化学求解器，而是在现有粒子追踪框架中，以尽量小的改动实现：

1. 粒子携带可反应浓度
2. 该浓度随粒子在裂缝中的停留时间而衰减
3. 浓度损失转换为反应物消耗量
4. 反应物消耗量通过化学计量和矿物摩尔体积转换为矿物体积变化
5. 矿物体积变化进一步转换为裂缝开度变化

因此，这一模块本质上是一个：

`粒子输运 + 一阶反应衰减 + 质量守恒几何更新`

的耦合模块。

## 3. 当前实现的核心思想

### 3.1 粒子携带反应浓度

在 [Particle.h](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/Transport/Particle.h) 中，每个粒子包含两个与 chemistry 直接相关的量：

- `reactive_concentration`
- `representative_volume`

含义分别是：

- `reactive_concentration`
  粒子当前携带的可反应浓度

- `representative_volume`
  该数值粒子所代表的真实流体体积

其中第二项非常关键，因为从浓度变化转换为摩尔数变化时，必须乘以体积。

### 3.2 粒子浓度的时间演化

在 [Chemistry.cpp](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/Chemistry/Chemistry.cpp) 中，浓度演化由：

`UpdateReactiveConcentration(double concentration_in, double residence_time)`

实现。

当前使用的一阶衰减关系为：

```text
C_out = C_in * exp(-k * dt)
```

其中：

- `C_in` 是进入当前段前的粒子浓度
- `C_out` 是离开当前段或本次反应步后的粒子浓度
- `k` 是 `REACTIVE_CONCENTRATION_DECAY`
- `dt` 是该粒子在当前裂缝段中的停留时间

这个关系的优点是：

- 数值稳定
- 不会自然产生负浓度
- 易于和粒子停留时间耦合

但它是一个简化反应律，并不包含：

- 饱和指数
- 平衡浓度
- pH
- 温度
- 表面积动力学
- 多组分耦合

因此，它应理解为“有效反应能力衰减模型”，而不是完整地球化学模型。

## 4. 质量守恒如何实现

### 4.1 从浓度损失到反应物摩尔数

在 `EvaluateReactiveStep(...)` 中，首先计算浓度损失：

```text
ΔC = C_in - C_out
```

然后乘以该粒子代表的流体体积：

```text
Δn_reactive = ΔC * V_particle
```

其中：

- `Δn_reactive` 的单位是 `mol`
- `V_particle` 即 `representative_volume`

这一步是质量守恒的关键桥梁。

### 4.2 从反应物摩尔数到矿物体积变化

随后代码使用：

```text
Δn_mineral = Δn_reactive * ν
ΔV_mineral = Δn_mineral * V_m
```

其中：

- `ν = REACTIVE_TO_MINERAL_STOICH`
- `V_m = MINERAL_MOLAR_VOLUME`

因此：

- `REACTIVE_TO_MINERAL_STOICH` 控制化学计量关系
- `MINERAL_MOLAR_VOLUME` 控制每摩尔矿物对应多少固体体积变化

### 4.3 从矿物体积变化到裂缝开度变化

当前模型是二维 DFN，因此需要一个出平面厚度：

- `FRACTURE_OUT_OF_PLANE_THICKNESS`

局部开度变化写为：

```text
Δb_local = ΔV_mineral / (L_traversed * thickness)
```

因为当前代码中每条 fracture segment 只存储一个统一 aperture，而不是子段 aperture 场，所以最终应用的是等效段平均变化：

```text
Δb_segment = (L_traversed / L_segment) * Δb_local
```

这等价于：

```text
Δb_segment = ΔV_mineral / (L_segment * thickness)
```

## 5. 为什么粒子代表体积是动态的

当前代码不是把每个粒子赋予固定体积，而是在注入时根据入口流动状态计算：

```text
V_particle = |u_inlet| * aperture_inlet * thickness * dt_particle
```

这在 [Transport.cpp](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/Transport/Transport.cpp) 中的

- `ComputeParticleRepresentativeVolumeAtInjection(...)`

实现。

它的物理含义是：

每个粒子代表“在该注入时间份额内，通过入口裂缝的流体体积”。

这样做比固定常数体积更适合压力差边界条件，因为入口流速和开度会随网络状态变化。

## 6. 局部更新与部分穿越

当前实现支持“一个 reaction/update interval 内，粒子只走过裂缝段的一部分”这种情况。

代码中会记录：

- `segment_length`
- `traversed_length`

如果粒子没有走完整段，则只对走过的那部分计算局部反应，并折算为该段的等效平均 aperture change。

因此，当前实现不是简单地“粒子碰到该段就给整段一次完整更新”，而是考虑了实际走过的距离比例。

## 7. 零浓度粒子的处理

如果粒子浓度降到零或非常接近零：

- 该粒子仍然可以继续输运
- 但它不再继续改几何

这在逻辑上表示：

- 它仍然是一个示踪/流体包
- 但它携带的可反应组分已经耗尽

这是当前代码中非常重要的设计选择。

## 8. 当前 chemistry 参数是如何进入程序的

现在 chemistry 参数可以通过：

- `Input/Chemistry_files/*.txt`

输入。

由 [Parameters.cpp](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/Input_Output/Parameters.cpp) 读取后，
在 [PERFORM.cpp](/home/zhouw/Documents/Codes/ReactiveSoluteTransport/Code/src/PERFORM.cpp) 中通过：

- `ConfigureChemistryParameters(...)`

写入运行时 chemistry 模块。

因此，当前程序已经支持：

- 不改代码
- 只换 chemistry 输入文件

来运行不同反应性 case。

## 9. 当前模型的物理背景与适用解释

从原理上看，当前化学模块适合解释为：

“粒子携带有限可反应物，沿裂缝输运时逐步消耗，该消耗通过质量守恒映射到矿物体积变化，并最终导致裂缝开度演化。”

它更接近以下类型的简化模型：

- 反应能力衰减模型
- 有效一阶动力学模型
- 质量守恒驱动的几何演化模型

而不属于：

- 完整 PHREEQC/PFLOTRAN 型多组分平衡-动力学反应模型
- 具有显式饱和度控制的矿物反应模型

## 10. 当前实现中“已实现”的部分

明确已实现：

1. 粒子携带浓度
2. 浓度按停留时间衰减
3. 段内局部反应更新
4. 部分穿越比例处理
5. 浓度损失到摩尔数的质量守恒转换
6. 摩尔变化到矿物体积变化的转换
7. 矿物体积变化到 aperture change 的转换
8. 更新后的 aperture 反馈到后续流场与输运

## 11. 当前实现中“未实现或仅简化处理”的部分

未实现或只做了简化：

1. 多组分反应网络
2. 显式平衡浓度/饱和指数
3. pH、温度、离子强度影响
4. 基于表面积或矿物反应面积的动力学
5. 真实矿物学相变追踪
6. 逐子段 aperture 场
   当前仍是“每个 fracture segment 一个 aperture”
7. 直接导出的粒子浓度历史
   当前浓度可由模型公式重构，但默认没有逐粒子浓度输出文件

## 12. 如何理解不同 chemistry case

当前 chemistry case 的差异，本质上体现在两个核心量：

- `C0`
- `k_decay`

因为当前小步近似下，几何响应强度大体与：

```text
C0 * k_decay
```

成比例。

因此：

- 增大 `C0`
  会增强每个粒子的反应物库存

- 增大 `k_decay`
  会加快反应物随停留时间的消耗

二者共同决定：

- 初始反应强度
- 反应持续时间
- 下游几何更新还能剩多少“化学驱动力”

## 13. 总结

当前代码中的 chemistry 模块，是一个面向 DFN 粒子输运框架设计的、模块化的一阶反应-几何耦合模块。

它的最大特点是：

- 保持了粒子法输运主框架
- 用质量守恒把浓度变化和几何变化连接起来
- 支持通过输入文件快速切换不同 chemistry case

如果后续希望它更接近真实地球化学模型，最自然的下一步扩展方向是：

1. 输出并分析逐粒子浓度历史
2. 引入平衡浓度 `C_eq`
3. 将反应速率改写为 `dC/dt = -k(C-C_eq)` 一类形式
4. 进一步加入更真实的矿物反应速率律
