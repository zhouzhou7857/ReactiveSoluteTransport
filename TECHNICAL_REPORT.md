# Implementation-Based Technical Report

## Scope

This report describes the current active code under `Code/src` and the maintained input format under `Input/`.

## A. Overall purpose

The code simulates particle-based solute transport in a two-dimensional discrete fracture network embedded in a porous matrix. The current main path couples transport to fracture evolution through a PHREEQC-calibrated cumulative mineral-volume law:

- each injected particle carries a representative fluid volume;
- each particle also stores a cumulative reaction age `particle_age`;
- the chemistry module evaluates the cumulative reference law `DeltaV_ref(t)`;
- the mineral-volume increment over one transport interval is
  `DeltaV_ref(t_end) - DeltaV_ref(t_start)`;
- that increment is scaled by `V_particle / Vref`, with optional additional geometry-based factors;
- the accumulated mineral volume on each segment is converted to a segment-averaged aperture increment;
- updated aperture modifies transmissivity and velocity;
- the flow field is recomputed when required.

The older concentration-decay path and empirical `Nb/Nt` aperture laws still exist as legacy code, but they are not the default mechanism on the current main execution path.

## B. Main workflow

### Entry and orchestration

Primary entry:

- `Code/src/PERFORM.cpp`

Main supporting modules:

- `Code/src/Input_Output/Parameters.h/.cpp`
- `Code/src/Domain_Definition/NetworkMeshes.h/.cpp`
- `Code/src/Domain_Definition/DFNComputation.h/.cpp`
- `Code/src/Transport/Transport.h/.cpp`
- `Code/src/Chemistry/Chemistry.h/.cpp`
- `Code/src/Input_Output/Results.h/.cpp`

### Active execution sequence

1. `PERFORM.cpp` constructs `Parameters param`.
2. `Parameters::Parameters()` reads `Input/File_names.txt`, or an override path passed on the command line.
3. `Parameters::read_param()` reads:
   - domain parameters;
   - simulation parameters;
   - DFN mode and DFN parameters.
4. `PERFORM.cpp` resets chemistry defaults, then configures the active chemistry law from:
   - simulation input (`Vref`, `fracture thickness`);
   - optional runtime environment overrides such as `RST_CHEM_MODE`, `RST_DELTA_V_A*`, `RST_DELTA_V_K*`, `RST_DELTA_V_L`, `RST_USE_EFFECTIVE_DIFFUSION_HEIGHT_FACTOR`, and `RST_USE_VP_WIDTH_CORRECTION`.
5. `PERFORM.cpp` selects the DFN construction path:
   - `generation_realistic3/4`:
     `NetworkMeshes(param.density_param, param.exponent_param, param, domain)`
   - all other modes:
     `NetworkMeshes(param.code_path, param.file_name_DFN, domain)`
6. The raw generated network is written to `Output/DFN_raw.txt`.
7. If `backbone=true`, `return_backbone()` extracts the connected flowing network.
8. `Transport::Particles_Transport()` performs stepwise transport and geometry update.
9. `Results` writes `cdf.txt` and `pdf.txt`.
10. Final DFN files are written:
   - `DFN_init.txt`
   - `DFN.txt`
   - `DFN_aperture_delta.txt`

## C. Input format

### File list

`Input/File_names.txt` is now used primarily as a 3-line file:

1. domain file
2. simulation file
3. DFN file

A fourth chemistry line is still accepted for backward compatibility, but the active main path does not read chemistry-input files anymore.

### Simulation file

The simulation file now includes two active chemistry-related inputs at the end:

11. `Vref` `[m^3]`
12. `fracture thickness` `[m]`

These are read in `Parameters.cpp` and then applied in `PERFORM.cpp`.

## D. DFN construction path

### Current path split

The key switch from the old `generation_realistic3` path to the current one is implemented by two core changes:

1. `Parameters.cpp` reads the parameter set required by `generation_realistic3/4` into `Parameters`.
2. `PERFORM.cpp` routes `generation_realistic3/4` to the constructor-based `NetworkMeshes(double,double,Parameters,Domain)` path instead of the older `NetworkMeshes(code_path,file_name_DFN,domain)` path.

### Current `generation_realistic3` behavior

The active implementation is in the constructor-based branch of `NetworkMeshes.cpp`.

For `generation_realistic3`:

- the loop stops when accumulated density reaches `density_param`;
- fracture center is uniformly sampled in the domain;
- fracture length is generated from `exponent_param`;
- angle is uniformly random in `[0, pi]`;
- aperture is sampled from a truncated lognormal distribution controlled by:
  - `b_min`
  - `b_max`
  - `mean_lnb`
  - `RSD_lnb`
- the current implementation uses the same random draw for both length and aperture sampling, so longer fractures tend to receive larger apertures.

Important note:

- `density_param` is not an exact fracture count target.
- the current `generation_realistic3` path is density-threshold-driven, not count-driven.

## E. Particle model

Each particle stores:

- current position;
- current mesh index;
- previous mesh index;
- mesh and intersection history;
- current time;
- injection time;
- cumulative distance in the current fracture;
- cumulative residence time in the current fracture;
- cumulative chemistry age `particle_age`;
- representative fluid volume `representative_volume`.

Injection assigns:

- inlet position and mesh;
- reset transport counters;
- representative volume based on inlet flux and particle time share:
  `V_particle = |u_inlet| * aperture_inlet * thickness * dt_particle`.

## F. Active chemistry-driven geometry update

### Active mechanism

The active geometry update is now based on the cumulative mineral-volume law in `Chemistry.cpp`.

The reference law is:

`DeltaV_ref(t) = A1 * (1 - exp(-k1 * t)) + A2 * (1 - exp(-k2 * t)) + L * t`

For one particle over one transport interval:

`DeltaV_particle[t_start,t_end] = (DeltaV_ref(t_end) - DeltaV_ref(t_start)) * scale`

where `scale` is:

- always `V_particle / Vref`;
- optionally multiplied by an effective diffusion-height factor;
- optionally multiplied by a Vp-width correction factor.

`Transport.cpp` accumulates these mineral-volume increments per segment, then applies the total segment change through:

`delta_b = delta_V / (segment_length * thickness * 2)`

### Optional scaling controls

The current code supports optional runtime controls:

- `RST_USE_EFFECTIVE_DIFFUSION_HEIGHT_FACTOR`
- `RST_EFFECTIVE_DIFFUSION_COEFFICIENT`
- `RST_EFFECTIVE_DIFFUSION_TIME`
- `RST_USE_VP_WIDTH_CORRECTION`

These are intended for benchmark testing and sensitivity analysis. They are not mandatory on the default path.

## G. Legacy code status

The following legacy interfaces still exist:

- `UpdateReactiveConcentration()`
- `EvaluateReactiveStep()`
- `ComputeAperture()`
- `ComputeDeltaAperture()`
- `NetworkMeshes::ChangeAperture()`
- `NetworkMeshes::ChangeApertureByRatio()`

They are retained for compatibility and historical reference, but they are not the default mechanism on the current main execution path.

## H. Outputs

The current run can produce:

- `Output/DFN_raw.txt`
- `Output/DFN_init.txt`
- `Output/DFN.txt`
- `Output/DFN_aperture_delta.txt`
- `Output/DFN_step*.txt`
- `Output/particle_positions_t*.csv`
- `Output/cdf.txt`
- `Output/pdf.txt`

`DFN_raw.txt` is especially useful for inspecting the generated DFN before backbone filtering.
