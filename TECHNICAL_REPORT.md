# Implementation-Based Technical Report

## Scope

This report reflects the active code under `Code/src` and the currently maintained input format under `Input/`.

## A. Overall purpose

The code simulates particle-based solute transport in a two-dimensional discrete fracture network embedded in a porous matrix. The active main path couples transport to DFN evolution through chemistry-driven aperture updates:

- each injected particle carries an initial reactive concentration;
- concentration decays during residence inside a fracture segment;
- reactive mass loss is converted to mineral volume change;
- mineral volume change is converted to a segment-averaged aperture increment;
- updated aperture modifies transmissivity and velocity;
- the flow field is recomputed when needed.

This is different from the older empirical particle-throughput aperture law, which is still present as legacy code but is not the default active mechanism in the current main path.

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
2. `Parameters::Parameters()` reads `Input/File_names.txt`.
3. `Parameters::read_param()` reads:
   - domain parameters;
   - simulation parameters;
   - DFN mode and DFN parameters;
   - optional chemistry parameters.
4. `PERFORM.cpp` configures active chemistry constants through `ConfigureChemistryParameters(...)`.
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

## C. DFN construction path

### Current path split

The key switch from the old `generation_realistic3` path to the current one is implemented by two core changes:

1. `Parameters.cpp` now reads the parameter set required by `generation_realistic3/4` into `Parameters`.
2. `PERFORM.cpp` now routes `generation_realistic3/4` to the constructor-based `NetworkMeshes(double,double,Parameters,Domain)` path instead of the older `NetworkMeshes(code_path,file_name_DFN,domain)` path.

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

## D. Particle model

Each particle stores:

- current position;
- current mesh index;
- previous mesh index;
- mesh and intersection history;
- current time;
- injection time;
- cumulative distance in the current fracture;
- cumulative residence time in the current fracture;
- reactive concentration carried by the particle;
- representative fluid volume carried by the particle.

Injection assigns:

- inlet position and mesh;
- reset transport counters;
- `reactive_concentration = INITIAL_REACTIVE_CONCENTRATION`;
- representative volume based on inlet flux and particle time share.

## E. Chemistry-driven geometry update

### Active mechanism

The active geometry update no longer relies on the old empirical `Nb/Nt` aperture law on the main path.

Instead:

1. `UpdateReactiveConcentration()` computes:
   - `C_out = C_in * exp(-k * t_residence)`
2. `EvaluateReactiveStep()` computes:
   - reactive moles consumed;
   - mineral moles changed;
   - mineral volume changed;
   - local aperture change;
   - segment-averaged aperture change.
3. `Transport.cpp` accumulates those changes in `step_aperture_change`.
4. Segment apertures are updated by adding the accumulated segment-averaged increment.
5. If aperture reaches zero, the network may be rebuilt and flow recomputed.

### Required chemistry parameters

The current chemistry file provides:

- initial reactive concentration;
- concentration decay rate;
- reactive-to-mineral stoichiometric factor;
- mineral molar volume;
- fracture out-of-plane thickness.

## F. Legacy code status

The following legacy interfaces still exist:

- `ComputeAperture()`
- `ComputeDeltaAperture()`
- `NetworkMeshes::ChangeAperture()`
- `NetworkMeshes::ChangeApertureByRatio()`

They are retained for compatibility and historical reference, but they are not the default mechanism on the current main execution path.

## G. Outputs

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
