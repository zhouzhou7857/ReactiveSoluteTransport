# Validation report

## Active chemistry implementation
- `ComputeParticleSegmentMineralVolumeChange` uses `delta_v_ref * scale`: `True`
- `scale` is computed from `V_particle / V_ref` when optional corrections are disabled.

## Volume checks
### t_cross = 1e+02 s
- `Q_in = |u| * b * thickness = 5.000000000000e-07 m^3/s`
- `V_ref = Q_in * t_injection = 5.000000000000e-02 m^3`
- `particle_time_share = 1.000100010001e+01 s`
- `nominal V_particle = 5.000500050005e-06 m^3`
- `sum(V_particle)_nominal = 5.000500050005e-02 m^3`
- `sum(V_particle)/V_ref = 1.000100010001e+00`
- `PHREEQC water = 5.000000000000e+01 kg`
- `PHREEQC A_ref = 2.000000000000e-02 m^2`
- `one F_ref(t) file shared across segment counts = True`
### t_cross = 1e+03 s
- `Q_in = |u| * b * thickness = 5.000000000000e-08 m^3/s`
- `V_ref = Q_in * t_injection = 5.000000000000e-03 m^3`
- `particle_time_share = 1.000100010001e+01 s`
- `nominal V_particle = 5.000500050005e-07 m^3`
- `sum(V_particle)_nominal = 5.000500050005e-03 m^3`
- `sum(V_particle)/V_ref = 1.000100010001e+00`
- `PHREEQC water = 5.000000000000e+00 kg`
- `PHREEQC A_ref = 2.000000000000e-02 m^2`
- `one F_ref(t) file shared across segment counts = True`
### t_cross = 1e+04 s
- `Q_in = |u| * b * thickness = 5.000000000000e-09 m^3/s`
- `V_ref = Q_in * t_injection = 5.000000000000e-04 m^3`
- `particle_time_share = 1.000100010001e+01 s`
- `nominal V_particle = 5.000500050005e-08 m^3`
- `sum(V_particle)_nominal = 5.000500050005e-04 m^3`
- `sum(V_particle)/V_ref = 1.000100010001e+00`
- `PHREEQC water = 5.000000000000e-01 kg`
- `PHREEQC A_ref = 2.000000000000e-02 m^2`
- `one F_ref(t) file shared across segment counts = True`
### t_cross = 1e+05 s
- `Q_in = |u| * b * thickness = 5.000000000000e-10 m^3/s`
- `V_ref = Q_in * t_injection = 5.000000000000e-05 m^3`
- `particle_time_share = 1.000100010001e+01 s`
- `nominal V_particle = 5.000500050005e-09 m^3`
- `sum(V_particle)_nominal = 5.000500050005e-05 m^3`
- `sum(V_particle)/V_ref = 1.000100010001e+00`
- `PHREEQC water = 5.000000000000e-02 kg`
- `PHREEQC A_ref = 2.000000000000e-02 m^2`
- `one F_ref(t) file shared across segment counts = True`

## Workflow checks
- No benchmark-local PHREEQC input uses per-particle `V_particle` as the water reference volume.
- One PHREEQC file is generated per crossing-time case, not per particle and not per segment count.
- `A_ref` is fixed to the fracture two-wall area and is not multiplied by `Np`.
