# Single-Fracture Pure Gypsum Benchmark

This folder contains the automation for a clean benchmark of the current particle-tracking chemistry model on one straight fracture in a `10 cm x 10 cm` square domain.

The benchmark is designed to answer one specific question:

- for which residence-time scale, equivalently `Re` and `Da`, does the current segment-based implementation reproduce pure gypsum dissolution accurately enough;
- how strongly does the dissolved gypsum volume depend on fracture segmentation (`n = 1, 5, 10, 50, 100`).

Setup choices:

- one centered horizontal fracture spanning the full domain length;
- aperture fixed initially at `1 mm`;
- chemistry reduced to pure gypsum with the active law
  `DeltaV(t) = A1 * (1 - exp(-k1 * t))`;
- four target crossing times: `100 s`, `1000 s`, `1e4 s`, `1e5 s`;
- a separate output tree under `Output/benchmark_gypsum_single_fracture_10cm/`.

Commands:

```bash
python3 Benchmark/gypsum_single_fracture_10cm/run_benchmark.py
python3 Benchmark/gypsum_single_fracture_10cm/analyze_benchmark.py
```
