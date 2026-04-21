# Single-Fracture Pure Gypsum Benchmark

This folder contains the original `1 mm` single-fracture gypsum benchmark workflow.

It is retained as the baseline automation for:

- a `10 cm x 10 cm` square domain;
- one centered horizontal fracture spanning the full domain length;
- segment-count sweeps `n = 1, 5, 10, 50, 100`;
- crossing-time sweeps `100 s`, `1000 s`, `1e4 s`, `1e5 s`.

Current status:

- this benchmark is a preserved baseline, not the newest sensitivity workflow;
- the more recent `5 mm` benchmark lives under `Benchmark/gypsum_single_fracture_10cm_ap5mm/`;
- the newer `1 mm / 0.5 mm` threshold tests live under
  `Benchmark/gypsum_single_fracture_10cm_effdiff_threshold_n100/`.

The current maintained chemistry model in the codebase is the cumulative `DeltaV` law in `Code/src/Chemistry/Chemistry.cpp`. This older benchmark folder should therefore be read as an archived benchmark setup tied to its own output snapshot under `Output/benchmark_gypsum_single_fracture_10cm/`.
