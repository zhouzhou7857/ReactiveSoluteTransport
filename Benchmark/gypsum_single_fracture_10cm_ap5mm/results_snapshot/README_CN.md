# 5 mm 裂缝宽度版本 Benchmark

这组 benchmark 与当前 10 cm 单裂缝石膏 benchmark 的唯一差别是：

- 初始裂缝宽度改为 `5 mm`，原测试组为 `1 mm`。

其他设定保持一致：

- 域大小 `10 cm x 10 cm`；
- 一条居中水平裂缝贯穿整个域；
- out-of-plane thickness 为 `10 cm`；
- chemistry 仍然只保留纯 gypsum 的 `fast2` 项；
- 每个 crossing time 都使用其对应 `Vp` 下重新拟合的 gypsum-only `fast2` 公式，并取 `Vref = Vp`；
- 这一版结果不再加入额外的有效扩散高度因子；
- `relative_error` 现在定义为最大 aperture 增量相对于 code-consistent intrinsic gypsum aperture 增量的误差，其中 `r = 4e-5 mol/m^2/s`、`Vm = 7.4e-5 m^3/mol`；
- 穿越时间为 `100 s`、`1000 s`、`1e4 s`、`1e5 s`；
- segment 数为 `n = 1, 5, 10, 50, 100`。

相较于 1 mm 版本，5 mm 裂缝会改变流量、Re 数以及总注入水量，因此溶解体积和 aperture 演化都会系统变化。

文件说明：

- `generated_inputs/`：本组 benchmark 实际使用的输入文件副本；
- `results_snapshot/`：本组 benchmark 的 summary 和图件副本；
- `Output/benchmark_gypsum_single_fracture_10cm_ap5mm/`：原始求解器输出和分析结果。
