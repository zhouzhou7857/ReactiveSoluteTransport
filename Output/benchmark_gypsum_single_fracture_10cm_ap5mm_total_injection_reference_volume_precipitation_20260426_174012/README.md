# results_total_injection_reference_volume_precipitation_20260426_174012

## What changed relative to the old workflow
- The previous benchmark represented gypsum dissolution.
- This rerun represents kinetic gypsum precipitation from supersaturated Ca-SO4 water.
- It uses one total-injection reference system per crossing-time case.
- `V_ref` is the total injected solution volume computed from the PT inlet condition and injection duration.
- `A_ref` is fixed to the physical two-wall fracture area `2 * 10 cm * 10 cm = 200 cm^2`.
- PHREEQC water volume is set to `V_ref`, not to per-particle `V_particle`.
- Supersaturated solution is `Ca=30 mmol/kgw`, `S(6)=30 mmol/kgw`.
- Precipitation kinetics use `A=0.02 m2`, `k1=1 mol/m2/s`, `k2=5e-07 mol/m2/s`.
- PT chemistry uses `dV_particle = [F_ref(t_end) - F_ref(t_start)] * (V_particle / V_ref)`.
- `F_ref(t)` is negative for precipitation, so the existing aperture update decreases aperture without changing the geometry formula.

## Folder contents
- `copied_original_inputs/`: copies of the original benchmark inputs and scripts used for comparison.
- `new_inputs/`: regenerated PT input files for this rerun.
- `phreeqc/`: regenerated PHREEQC inputs, raw outputs, extracted curves, and fitted coefficients.
- `pt_cases/`: PT logs and raw outputs for all segment-count cases.
- `figures/`: final plots for the rerun and old-vs-new comparison.
- `figures/crossing_time_mechanism_n100.csv`: data behind the crossing-time mechanism plot.
- `validation/`: validation summary and formula checks.
- `scripts_snapshot/`: copies of the active benchmark script and chemistry-relevant source files.

## Active chemistry law
- PHREEQC simulates kinetic gypsum precipitation from supersaturated water.
- PHREEQC empirical rate is `rate = -A * (k1 * (sqrt(SR)-1)^10 + k2 * (sqrt(SR)-1)^2)` when `SR > 1`.
- The PHREEQC cumulative precipitated volume is fitted by the rational form:
- `F_ref(t)=Pinf*[1-w*(1+a10*t)^(-1/9)-(1-w)*(1+a2*t)^(-1)]`.
- `Pinf` is stored as a signed mineral-volume change. It is negative here because precipitation reduces aperture.
- PT applies the fitted law as `dV_particle = [F_ref(t_end)-F_ref(t_start)] * (V_particle/V_ref)`.

## Figure guide
- `fig_f_ref_curves.png`: signed PHREEQC reference curves for each crossing-time case.
- `fig_f_ref_fit_check.png`: PHREEQC data versus rational fitted `F_ref(t)`.
- `fig_aperture_profile_n100.png`: final aperture profile along the 10 cm fracture for `n=100`.
- `fig_aperture_evolution_n100.png`: mean and minimum aperture evolution for `n=100`.
- `fig_total_precipitated_volume_n100.png`: accumulated precipitated gypsum volume over simulation time.
- `fig_crossing_time_mechanism_n100.png`: explains the non-monotonic maximum aperture reduction by comparing `V_inj`, local precipitation completion, and their product.

## Crossing-time mechanism
- `V_inj = Q_in * t_injection` decreases as `t_cross` increases because the inlet velocity decreases.
- `eta_seg = P(t_cross/100) / Pinf` is the inlet-segment precipitation completion for the `n=100` aperture profile.
- `V_inj * eta_seg` is a local precipitation proxy. It peaks at `t_cross = 1e4 s`, matching the largest local aperture reduction.

## Reference-system summary
### t_cross = 1e+02 s
- `u_inlet = 1.000000000000e-03 m/s`
- `Q_in = 5.000000000000e-07 m^3/s`
- `t_injection = 1.000000000000e+05 s`
- `V_ref = 5.000000000000e-02 m^3`
- `A_ref = 2.000000000000e-02 m^2 = 200.000000 cm^2`
- `Maximum Ca/S-limited gypsum precipitate = 1.108500000000e-04 m^3`
- `Maximum Ca/S-limited gypsum moles = 1.500000000000e+00 mol`
- `Fit: Pinf = -5.958490314590e-05 m^3, w = 9.829208153779e-01, a10 = 4.842430557079e-04 s^-1, a2 = 1.916340855748e-04 s^-1`
### t_cross = 1e+03 s
- `u_inlet = 1.000000000000e-04 m/s`
- `Q_in = 5.000000000000e-08 m^3/s`
- `t_injection = 1.000000000000e+05 s`
- `V_ref = 5.000000000000e-03 m^3`
- `A_ref = 2.000000000000e-02 m^2 = 200.000000 cm^2`
- `Maximum Ca/S-limited gypsum precipitate = 1.108500000000e-05 m^3`
- `Maximum Ca/S-limited gypsum moles = 1.500000000000e-01 mol`
- `Fit: Pinf = -5.874782752693e-06 m^3, w = 9.770149403280e-01, a10 = 4.939306432555e-03 s^-1, a2 = 1.461607621155e-03 s^-1`
### t_cross = 1e+04 s
- `u_inlet = 1.000000000000e-05 m/s`
- `Q_in = 5.000000000000e-09 m^3/s`
- `t_injection = 1.000000000000e+05 s`
- `V_ref = 5.000000000000e-04 m^3`
- `A_ref = 2.000000000000e-02 m^2 = 200.000000 cm^2`
- `Maximum Ca/S-limited gypsum precipitate = 1.108500000000e-06 m^3`
- `Maximum Ca/S-limited gypsum moles = 1.500000000000e-02 mol`
- `Fit: Pinf = -5.827140713191e-07 m^3, w = 9.729585131695e-01, a10 = 5.022456863058e-02 s^-1, a2 = 1.196697670354e-02 s^-1`
### t_cross = 1e+05 s
- `u_inlet = 1.000000000000e-06 m/s`
- `Q_in = 5.000000000000e-10 m^3/s`
- `t_injection = 1.000000000000e+05 s`
- `V_ref = 5.000000000000e-05 m^3`
- `A_ref = 2.000000000000e-02 m^2 = 200.000000 cm^2`
- `Maximum Ca/S-limited gypsum precipitate = 1.108500000000e-07 m^3`
- `Maximum Ca/S-limited gypsum moles = 1.500000000000e-03 mol`
- `Fit: Pinf = -6.119775374372e-08 m^3, w = 9.737179551403e-01, a10 = 5.254735011434e-01 s^-1, a2 = 1.669979509308e-27 s^-1`

## n=100 final results
### t_cross = 1e+02 s
- `mineral_volume_change = -3.253569400000e-07 m^3`
- `precipitated_volume = 3.253569400000e-07 m^3`
- `min_delta_aperture = -1.670040000000e-02 mm`
- `mean_delta_aperture = -1.626784700000e-02 mm`
### t_cross = 1e+03 s
- `mineral_volume_change = -1.094278720000e-06 m^3`
- `precipitated_volume = 1.094278720000e-06 m^3`
- `min_delta_aperture = -1.604360000000e-01 mm`
- `mean_delta_aperture = -5.471393600000e-02 mm`
### t_cross = 1e+04 s
- `mineral_volume_change = -2.939198140000e-07 m^3`
- `precipitated_volume = 2.939198140000e-07 m^3`
- `min_delta_aperture = -5.409960000000e-01 mm`
- `mean_delta_aperture = -1.469599070000e-02 mm`
### t_cross = 1e+05 s
- `mineral_volume_change = -3.962873114700e-08 m^3`
- `precipitated_volume = 3.962873114700e-08 m^3`
- `min_delta_aperture = -1.491250000000e-01 mm`
- `mean_delta_aperture = -1.981436557350e-03 mm`
