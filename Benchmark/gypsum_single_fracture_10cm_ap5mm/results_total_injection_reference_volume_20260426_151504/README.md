# results_total_injection_reference_volume_20260426_151504

## What changed relative to the old workflow
- The old benchmark used one particle-volume-specific cumulative law `F_Vp(t)` per crossing-time case and set `Vref = Vp`.
- This rerun uses one total-injection reference system per crossing-time case.
- `V_ref` is the total injected solution volume computed from the PT inlet condition and injection duration.
- `A_ref` is fixed to the physical two-wall fracture area `2 * 10 cm * 10 cm = 200 cm^2`.
- PHREEQC water volume is set to `V_ref`, not to per-particle `V_particle`.
- PT chemistry uses `dV_particle = [F_ref(t_end) - F_ref(t_start)] * (V_particle / V_ref)`.

## Folder contents
- `copied_original_inputs/`: copies of the original benchmark inputs and scripts used for comparison.
- `new_inputs/`: regenerated PT input files for this rerun.
- `phreeqc/`: regenerated PHREEQC inputs, raw outputs, extracted curves, and fitted coefficients.
- `pt_cases/`: PT logs and raw outputs for all segment-count cases.
- `figures/`: final plots for the rerun and old-vs-new comparison.
- `validation/`: validation summary and formula checks.
- `scripts_snapshot/`: copies of the active benchmark script and chemistry-relevant source files.

## Reference-system summary
### t_cross = 1e+02 s
- `u_inlet = 1.000000000000e-03 m/s`
- `Q_in = 5.000000000000e-07 m^3/s`
- `t_injection = 1.000000000000e+05 s`
- `V_ref = 5.000000000000e-02 m^3`
- `A_ref = 2.000000000000e-02 m^2 = 200.000000 cm^2`
- `Gypsum volume = 5.000000000000e-02 m^3`
- `Gypsum moles = 6.765899864682e+02 mol`
- `Fit: A2 = 1.685800000000e-05 m^3, k2 = 1.000000000000e-03 s^-1, L = 1.000000000000e-10 m^3/s`
### t_cross = 1e+03 s
- `u_inlet = 1.000000000000e-04 m/s`
- `Q_in = 5.000000000000e-08 m^3/s`
- `t_injection = 1.000000000000e+05 s`
- `V_ref = 5.000000000000e-03 m^3`
- `A_ref = 2.000000000000e-02 m^2 = 200.000000 cm^2`
- `Gypsum volume = 5.000000000000e-03 m^3`
- `Gypsum moles = 6.765899864682e+01 mol`
- `Fit: A2 = 5.470700000000e-06 m^3, k2 = 1.000000000000e-03 s^-1, L = 1.000000000000e-10 m^3/s`
### t_cross = 1e+04 s
- `u_inlet = 1.000000000000e-05 m/s`
- `Q_in = 5.000000000000e-09 m^3/s`
- `t_injection = 1.000000000000e+05 s`
- `V_ref = 5.000000000000e-04 m^3`
- `A_ref = 2.000000000000e-02 m^2 = 200.000000 cm^2`
- `Gypsum volume = 5.000000000000e-04 m^3`
- `Gypsum moles = 6.765899864682e+00 mol`
- `Fit: A2 = 5.562200000000e-07 m^3, k2 = 1.000000000000e-04 s^-1, L = 1.000000000000e-10 m^3/s`
### t_cross = 1e+05 s
- `u_inlet = 1.000000000000e-06 m/s`
- `Q_in = 5.000000000000e-10 m^3/s`
- `t_injection = 1.000000000000e+05 s`
- `V_ref = 5.000000000000e-05 m^3`
- `A_ref = 2.000000000000e-02 m^2 = 200.000000 cm^2`
- `Gypsum volume = 5.000000000000e-05 m^3`
- `Gypsum moles = 6.765899864682e-01 mol`
- `Fit: A2 = 5.562300000000e-08 m^3, k2 = 1.000000000000e-05 s^-1, L = 1.000000000000e-10 m^3/s`

## n=100 final results
### t_cross = 1e+02 s
- `dissolved_volume = 1.627906800000e-06 m^3`
- `max_delta_aperture = 8.547940000000e-02 mm`
- `mean_delta_aperture = 8.139534000000e-02 mm`
### t_cross = 1e+03 s
- `dissolved_volume = 3.658923200000e-06 m^3`
- `max_delta_aperture = 2.850630000000e-01 mm`
- `mean_delta_aperture = 1.829461600000e-01 mm`
### t_cross = 1e+04 s
- `dissolved_volume = 1.296932760000e-06 m^3`
- `max_delta_aperture = 7.828760000000e-02 mm`
- `mean_delta_aperture = 6.484663800000e-02 mm`
### t_cross = 1e+05 s
- `dissolved_volume = 5.186120260000e-06 m^3`
- `max_delta_aperture = 5.286860000000e-01 mm`
- `mean_delta_aperture = 2.593060130000e-01 mm`
