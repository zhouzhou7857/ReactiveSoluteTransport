# Single-Fracture Pure Gypsum Benchmark: 5 mm Aperture Variant

This benchmark is identical to the current 10 cm single-fracture gypsum test except for one change:

- initial fracture aperture is `5 mm` instead of `1 mm`.

Shared setup:

- domain size `10 cm x 10 cm`;
- one centered horizontal fracture spanning the full domain length;
- out-of-plane thickness `10 cm`;
- chemistry kept as pure gypsum using only the `fast2` term;
- each crossing-time case uses its own gypsum-only `fast2` law fitted for the corresponding particle volume `Vp`, with `Vref = Vp`;
- no additional effective diffusion-height factor is applied in this result set;
- `relative_error` is defined from the maximum simulated aperture increase relative to the code-consistent intrinsic gypsum aperture increase using `r = 4e-5 mol/m^2/s` and `Vm = 7.4e-5 m^3/mol`;
- crossing times `100 s`, `1000 s`, `1e4 s`, `1e5 s`;
- segment counts `n = 1, 5, 10, 50, 100`.

Compared with the 1 mm benchmark, the larger aperture changes both flow rate and Reynolds number.

Files:

- `generated_inputs/`: copies of the exact input files used for this variant.
- `results_snapshot/`: copied summaries and figures for this benchmark.
- `Output/benchmark_gypsum_single_fracture_10cm_ap5mm/`: raw solver outputs and analysis products.
