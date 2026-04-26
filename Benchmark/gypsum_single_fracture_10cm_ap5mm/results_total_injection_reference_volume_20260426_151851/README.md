# results_total_injection_reference_volume_20260426_151851

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
- `Fit: A2 = 2.131932783233e-04 m^3, k2 = 2.774203696879e-07 s^-1, L = 2.711015461813e-19 m^3/s`
### t_cross = 1e+03 s
- `u_inlet = 1.000000000000e-04 m/s`
- `Q_in = 5.000000000000e-08 m^3/s`
- `t_injection = 1.000000000000e+05 s`
- `V_ref = 5.000000000000e-03 m^3`
- `A_ref = 2.000000000000e-02 m^2 = 200.000000 cm^2`
- `Gypsum volume = 5.000000000000e-03 m^3`
- `Gypsum moles = 6.765899864682e+01 mol`
- `Fit: A2 = 7.270810710158e-06 m^3, k2 = 8.193391570945e-06 s^-1, L = 8.396448112282e-27 m^3/s`
### t_cross = 1e+04 s
- `u_inlet = 1.000000000000e-05 m/s`
- `Q_in = 5.000000000000e-09 m^3/s`
- `t_injection = 1.000000000000e+05 s`
- `V_ref = 5.000000000000e-04 m^3`
- `A_ref = 2.000000000000e-02 m^2 = 200.000000 cm^2`
- `Gypsum volume = 5.000000000000e-04 m^3`
- `Gypsum moles = 6.765899864682e+00 mol`
- `Fit: A2 = 5.722657481985e-07 m^3, k2 = 1.053529781365e-04 s^-1, L = 1.982637032109e-31 m^3/s`
### t_cross = 1e+05 s
- `u_inlet = 1.000000000000e-06 m/s`
- `Q_in = 5.000000000000e-10 m^3/s`
- `t_injection = 1.000000000000e+05 s`
- `V_ref = 5.000000000000e-05 m^3`
- `A_ref = 2.000000000000e-02 m^2 = 200.000000 cm^2`
- `Gypsum volume = 5.000000000000e-05 m^3`
- `Gypsum moles = 6.765899864682e-01 mol`
- `Fit: A2 = 5.633450071382e-08 m^3, k2 = 1.079500437139e-03 s^-1, L = 1.102025953896e-39 m^3/s`

## n=100 final results
### t_cross = 1e+02 s
- `dissolved_volume = 5.911905800000e-09 m^3`
- `max_delta_aperture = 2.957590000000e-04 mm`
- `mean_delta_aperture = 2.955952900000e-04 mm`
### t_cross = 1e+03 s
- `dissolved_volume = 5.905925000000e-08 m^3`
- `max_delta_aperture = 2.979690000000e-03 mm`
- `mean_delta_aperture = 2.952962500000e-03 mm`
### t_cross = 1e+04 s
- `dissolved_volume = 3.586201180000e-07 m^3`
- `max_delta_aperture = 3.008010000000e-02 mm`
- `mean_delta_aperture = 1.793100590000e-02 mm`
### t_cross = 1e+05 s
- `dissolved_volume = 5.609476934449e-08 m^3`
- `max_delta_aperture = 1.879110000000e-01 mm`
- `mean_delta_aperture = 2.804738467225e-03 mm`
