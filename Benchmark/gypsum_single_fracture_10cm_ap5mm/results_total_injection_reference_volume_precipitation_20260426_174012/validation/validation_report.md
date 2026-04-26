# Validation report

## Active chemistry implementation
- `ComputeParticleSegmentMineralVolumeChange` uses `delta_v_ref * scale`: `True`
- `scale` is computed from `V_particle / V_ref` when optional corrections are disabled.
- Active benchmark chemistry mode is `rational_precipitation`.
- Active reference law is `F_ref(t)=Pinf*[1-w*(1+a10*t)^(-1/9)-(1-w)*(1+a2*t)^(-1)]`.
- `Pinf < 0` means precipitation; negative `F_ref` decreases aperture through the existing geometry update.

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
- `rational Pinf = -5.958490314590e-05 m^3`
- `rational w = 9.829208153779e-01`
- `rational a10 = 4.842430557079e-04 s^-1`
- `rational a2 = 1.916340855748e-04 s^-1`
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
- `rational Pinf = -5.874782752693e-06 m^3`
- `rational w = 9.770149403280e-01`
- `rational a10 = 4.939306432555e-03 s^-1`
- `rational a2 = 1.461607621155e-03 s^-1`
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
- `rational Pinf = -5.827140713191e-07 m^3`
- `rational w = 9.729585131695e-01`
- `rational a10 = 5.022456863058e-02 s^-1`
- `rational a2 = 1.196697670354e-02 s^-1`
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
- `rational Pinf = -6.119775374372e-08 m^3`
- `rational w = 9.737179551403e-01`
- `rational a10 = 5.254735011434e-01 s^-1`
- `rational a2 = 1.669979509308e-27 s^-1`
- `one F_ref(t) file shared across segment counts = True`

## Workflow checks
- No benchmark-local PHREEQC input uses per-particle `V_particle` as the water reference volume.
- One PHREEQC file is generated per crossing-time case, not per particle and not per segment count.
- `A_ref` is fixed to the fracture two-wall area and is not multiplied by `Np`.
- Positive PHREEQC precipitated gypsum volume is exported to PT as negative `F_ref(t)` so aperture decreases.
- The optional `V_p` width correction and effective diffusion height factor are disabled for this benchmark.
