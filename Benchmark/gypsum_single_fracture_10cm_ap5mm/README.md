# Single-Fracture Pure Gypsum Benchmark: 5 mm Aperture Variant

This benchmark is identical to the current 10 cm single-fracture gypsum test except for one change:

- initial fracture aperture is `5 mm` instead of `1 mm`.

Shared setup:

- domain size `10 cm x 10 cm`;
- one centered horizontal fracture spanning the full domain length;
- out-of-plane thickness `10 cm`;
- chemistry kept as pure gypsum using only the `fast2` term;
- all crossing-time cases share one gypsum-only volume law with `A2 = 2.50122e-6`, `k2 = 2.72730e-5`, and `Vref = 1e-3 m^3`;
- an effective diffusion-height factor `m = sqrt(D * t_cross) / b` is applied with `D = 1e-9 m^2/s` and `b` approximated by the initial aperture in the analytical reference plots;
- crossing times `100 s`, `1000 s`, `1e4 s`, `1e5 s`;
- segment counts `n = 1, 5, 10, 50, 100`.

Compared with the 1 mm benchmark, the larger aperture changes both flow rate and Reynolds number.

Files:

- `generated_inputs/`: copies of the exact input files used for this variant.
- `results_snapshot/`: copied summaries and figures for this benchmark.
- `Output/benchmark_gypsum_single_fracture_10cm_ap5mm/`: raw solver outputs and analysis products.
