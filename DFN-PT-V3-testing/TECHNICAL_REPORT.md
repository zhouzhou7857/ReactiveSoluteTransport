# Implementation-Based Technical Report

## Scope

This report is based strictly on the active code under `Code/src`. Backup directories, historical snapshots, and output-side analysis scripts are excluded except where they help identify what is no longer on the main execution path.

## A. Overall purpose of the code

The code implements particle-based solute transport in a two-dimensional discrete fracture network (DFN) embedded in a porous matrix. The active workflow couples transport to DFN evolution by updating fracture aperture from particle-derived passage metrics and then recomputing the hydraulic flow field on the modified network.

The implemented model is not a generic geochemical simulator. In the current code, the structural evolution is driven by an empirical aperture law linked to particle traffic through each fracture mesh. The main coupling therefore is:

- transport produces per-mesh particle-passage measures;
- those measures update fracture aperture;
- updated aperture changes hydraulic conductance and transmissivity;
- the DFN flow field is recomputed;
- the next transport step uses the updated DFN.

## B. Main workflow of the program

### Main entry files and modules

Primary entry and orchestration:

- `Code/src/PERFORM.cpp`

Core modules:

- `Code/src/Input_Output/Parameters.h/.cpp`
- `Code/src/Domain_Definition/Domain.h/.cpp`
- `Code/src/Domain_Definition/NetworkMeshes.h/.cpp`
- `Code/src/Domain_Definition/DFNComputation.h/.cpp`
- `Code/src/Domain_Definition/FractureMesh.h/.cpp`
- `Code/src/Transport/Transport.h/.cpp`
- `Code/src/Transport/Projection.cpp`
- `Code/src/Transport/Transfer.cpp`
- `Code/src/Chemistry/Chemistry.h/.cpp`
- `Code/src/Input_Output/Results.h/.cpp`

Auxiliary but not central to the active solver:

- `Code/src/Visualisation/*`
- `Code/src/RngStream/*`
- `Code/src/Utilitaries/*`

### Execution sequence

1. `main()` in `PERFORM.cpp` constructs `Parameters param` and `Domain domain(param.Lx, param.Ly)`.
2. `Parameters::Parameters()` reads `Input/File_names.txt`, then calls `read_param()`.
3. `Parameters::read_param()` reads:
   - domain size, matrix diffusion coefficient, porosity, and left/right heads;
   - particle count, transfer control, matrix option, BTC time range, random seed, injection duration, output interval, and reaction time step;
   - the DFN mode from the selected DFN file.
4. `NetworkMeshes init_net_mesh(param.code_path, param.file_name_DFN, domain)` reads or generates the DFN.
5. `return_backbone()` is called from `PERFORM.cpp` if `backbone=true`.
   - It removes disconnected components and low-velocity fractures.
   - It internally calls `ComputeFlowVelocities()`.
6. `Transport transport(rng_tracker, net_mesh, domain, param)` is created.
7. `transport.Particles_Transport(arrival_times, option_injection)` performs the active transport-evolution loop.
8. If transport succeeds:
   - `Results results(arrival_times, param);`
   - `results.post_processing();`
   - `results.writing(param.code_path);`
   - `transport.WritePositionSnapshotsCSV(param.code_path);`
9. Final DFN files are written:
   - current DFN,
   - aperture-delta DFN,
   - initial DFN.

### Active versus inactive transport paths

The active path is `Particles_Transport()` plus `particle_displacement_step()`.

Older routines such as `particle_displacement()` still exist and also contain aperture-update logic, but they are not used by the current main program path in `PERFORM.cpp`.

## C. Particle-tracking algorithm: implementation steps

### 1. Particle initialization

Particle state is defined in `Particle.h`. Each particle stores:

- current position `M`;
- current fracture mesh index `mesh_index`;
- previous mesh index `prev_mesh_index`;
- visited mesh history;
- visited intersection history;
- current time `t`;
- injection time `t_injection`;
- cumulative distance traveled in the current fracture `L_in_fract`;
- cumulative time spent in the current fracture `t_in_fract`.

Initialization is done in `Transport::Particles_Injection()`.

Implemented behavior:

- A vector of `nb_part` particles is created.
- All particles are initially unset spatially: `M = NOT_DEFINED`, `mesh_index = -1`.
- If `t_injection > 0`, particles are assigned staggered injection times over `[0, t_injection]`.
- If `t_injection <= 0`, all particles start at time `0`.

This means initialization is temporal first, spatial second. Actual inlet placement occurs later inside `Particles_Transport()` when a particle becomes ready to inject.

### 2. Particle injection

Injection logic is implemented in `Particles_Transport()`.

The code first collects active inlet positions from the current DFN by `CollectInputPositions()`. A mesh qualifies if one endpoint is on the left boundary.

Three injection modes are implemented:

- `option_injection = 0`: extended uniform injection over all inlet meshes.
- `option_injection = 1`: localized injection at the inlet mesh closest to the midpoint of the left boundary.
- `option_injection = 2`: extended flux-weighted injection, with weights proportional to `aperture * |velocity|` on inlet meshes.

Particle counts per inlet mesh are computed by `ComputeInjectionCounts()`, which:

- normalizes weights,
- allocates the integer floor of expected counts,
- distributes remaining particles by largest fractional remainder.

When a particle becomes eligible for release (`mesh_index == -1` and `t <= t_current`), it is assigned:

- an inlet position,
- an inlet mesh,
- reset local transport counters.

### 3. Particle motion inside a fracture

The active stepwise mover is `particle_displacement_step()`.

It advances one particle over one reaction-time target `t_target = pa.t + dt_step`.

Two transport modes are implemented:

- infinite matrix: `infinite_matrix_displacement_step()`
- finite matrix: `finite_matrix_displacement_step()`

#### Infinite-matrix mode

Implemented logic:

1. Find the downstream endpoint by `distance_from_M_to_extremity()`.
2. Compute advection travel time to that endpoint by `Mesh_Time_Advection()`.
3. Convert advection time to total residence time by `Get_Total_Time_From_Advec_Time1D()`.
4. If the full event would overshoot `t_target`, move only partway along the fracture using `Mesh_Scale_Advection()`.
5. Otherwise move to the fracture endpoint.

The implemented residence-time formula in `FractureMesh::Get_Total_Time_From_Advec_Time1D()` is:

- `B = sqrt(Dm) * porosity * advection_time / (0.5 * aperture)`
- `res_time = advection_time + 0.25 * (B / erfc_inv(u))^2`

where `u` is a uniform random number.

This is explicit in the code. The function name still says `1D`, so the active residence-time model is a 1D matrix-diffusion correction to advection time.

#### Finite-matrix mode

Implemented logic:

1. Determine the downstream endpoint `M_out`.
2. Construct orthogonal projections from the current fracture toward nearby fractures with `Orthogonal_Projection_M()`.
3. Compute a stable sub-step advection time by `Advection_Time_Computation()`, using the local transfer-probability limit `proba_transfer`.
4. Compute total residence time from advection time by `Get_Total_Time_From_Advec_Time1D()`.
5. Interpret the extra part as diffusive time `t_diff = total_time - t_advec`.
6. Test whether a transfer to a neighboring fracture occurs before `t_diff` expires using `Transfer_Probability_and_Transfer_Time()`.
7. If transfer occurs:
   - choose transfer time,
   - choose which barrier/fracture is hit,
   - project the particle to the new fracture location.
8. Otherwise keep the particle in the same fracture and continue.
9. If the full event would exceed `t_target`, truncate movement and stop exactly at `t_target`.

This is a genuine transport-with-transfer implementation, not just intersection-based switching.

### 4. Path selection at fracture intersections

When a particle reaches a mesh endpoint and has not exited the domain, the next fracture is chosen by `Mesh_Neighbour_Choice_And_Test()`.

Implemented procedure:

1. Determine which endpoint of the current mesh the particle is closest to.
2. Gather neighboring meshes attached to that node.
3. Compute a signed outgoing flux proxy for each candidate:
   - `velocity * aperture` if flow leaves the node through the candidate orientation,
   - `-velocity * aperture` if orientation is reversed.
4. Convert those raw values to probabilities by `select_neighbour_slr()`.

`select_neighbour_slr()` implements a streamline-routing rule:

- zero valid downstream paths:
  - if already at the right boundary, allow exit;
  - otherwise warn and fail;
- one path:
  - probability 1;
- two paths:
  - prefer straight continuation if it dominates,
  - otherwise split or choose the bifurcation depending on relative flow;
- three or more paths:
  - probabilities proportional to outgoing flux.

This is explicitly implemented. It is not a purely geometric shortest-path rule.

### 5. Intersection handling

There are two distinct mechanisms:

- topological intersection handling at fracture nodes via `Mesh_Neighbour_Choice_And_Test()`;
- matrix-mediated transfer handling in finite-matrix mode via orthogonal projections and first-passage distributions.

Therefore the code supports both:

- transport along connected fracture graph edges;
- off-graph transfer to nearby fractures through the matrix.

### 6. Travel time and residence time calculation

Implemented time components:

- advection time along a fracture segment: `Mesh_Time_Advection()`;
- total residence time with matrix-diffusion correction: `Get_Total_Time_From_Advec_Time1D()`;
- finite-matrix transfer probability/time:
  - `Transfer_Probability_Computation()`
  - `Transfer_Time_Computation()`
  - `Transfer_Time_FPTD()`
  - `Transfer_Time_Feller()`

The finite-matrix transfer module uses:

- a first-passage formulation when only one barrier is present;
- a Feller-based two-barrier formulation when two barriers are present.

This is explicit from `Transfer.cpp`.

### 7. Boundary handling

Implemented boundary logic:

- left boundary:
  - inlet detection and particle injection;
- right boundary:
  - particle exit condition via `domain.on_output_limit()`;
- top and bottom boundaries:
  - hydraulic boundary conditions are set to zero Neumann by default in the active flow path;
  - transport-side reflection routines exist in `Projection.cpp`, notably `ReflectionOnBorder()`.

However, the active finite-matrix stepping path shown in `finite_matrix_displacement_step()` does not clearly call `ReflectionOnBorder()`. Therefore reflective transport at top/bottom boundaries appears supported by helper routines, but it is not clearly part of the currently active main stepping path.

### 8. Particle exit and output recording

Exit handling is in `Particles_Transport()`:

- if the particle reaches the right boundary, its arrival time is stored in `arrival_times[pa.no]`;
- the particle is then deactivated by setting `pa.t = -1`.

Output recording:

- `RecordPositions(time, particles)` stores snapshots in memory;
- `WriteDFNSnapshot()` writes `DFN_step*.txt`;
- `WritePositionSnapshotsCSV()` writes `particle_positions_t*.csv`;
- `Results::post_processing()` computes BTC CDF/PDF;
- `Results::writing()` writes `cdf.txt` and `pdf.txt`.

## D. Key particle-related variables and their physical meanings

- `nb_part`
  - total number of particles injected in the simulation.
- `proba_transfer`
  - upper control for transfer probability in finite-matrix advection substepping.
- `simu_option`
  - transport mode selector:
    - `0`: infinite matrix
    - `1`: finite matrix
- `t_injection`
  - duration over which particles are injected.
- `output_interval`
  - interval for writing particle position and DFN snapshots.
- `reaction_dt`
  - fixed coupling interval for DFN updates and flow recomputation.
- `pa.t`
  - current particle time.
- `pa.t_injection`
  - injection time of particle `pa`.
- `pa.t_in_fract`
  - cumulative time spent in the current fracture since entry.
- `pa.L_in_fract`
  - cumulative advective distance traveled in the current fracture since entry.
- `moved_distances[mesh_id]`
  - total advective distance accumulated by a particle in a mesh during a reaction step.
- `step_full_count[mesh_id]`
  - integer number of full mesh-length traversals in the current reaction step.
- `step_partial_sum[mesh_id]`
  - accumulated fractional traversals in the current reaction step.
- `Nb_eff`
  - effective particle-passage measure for a mesh in one reaction step:
    - full traversals plus fractional traversals.
- `reference_injected_count_dt`
  - reference injected particle count over one reaction step, used to normalize aperture updates.

## E. How particle transport influences DFN structure

### What is explicitly implemented

#### Fracture aperture

Implemented: yes.

Where:

- `Transport::Particles_Transport()` computes `Nb_eff` for each mesh after each reaction step.
- `NetworkMeshes::ChangeAperture()` updates the mesh aperture.
- `Chemistry::ComputeAperture()` computes the new aperture.

Update law:

- `ratio = min(1, Nb / Nt)`
- `delta = K_REACTION * dt * ratio`
- `b_new = b_old - delta`

with `K_REACTION = -1.5e-10`.

Because `K_REACTION` is negative, the implemented formula increases aperture when `ratio > 0`. The comments mention gypsum dissolution, which is consistent with aperture opening rather than closure.

The update is:

- local: mesh by mesh;
- periodic in time: once per `reaction_dt`;
- cumulative over the run;
- threshold-sensitive because aperture is clamped to nonnegative and later tested for zero.

#### Porosity

Implemented: no dynamic porosity update.

Evidence:

- `porosity` is read from input and stored in `Parameters` / `PhysicsParam`.
- It is used in residence-time and transfer calculations.
- No code updates `porosity` after initialization.

#### Permeability

Implemented explicitly as a mutable field: no.

There is no separate permeability state variable. What exists is hydraulic conductance/transmissivity derived from aperture.

#### Transmissivity

Implemented explicitly as a persistent evolving field: no.

Implemented implicitly as a function of aperture: yes.

Evidence:

- `HydraulicProperties::ReturnTransmissivity()` uses the cubic law `~ aperture^3`.
- The flow system is rebuilt from the current aperture.
- Therefore transmissivity changes indirectly after aperture changes, but the code does not store and update a separate transmissivity field mesh by mesh over time.

#### Local geometry

Implemented as coordinate movement or deformation: no.

Fracture endpoint coordinates are not updated by reaction. Aperture changes, but the segment geometry itself is fixed.

#### Fracture segment properties

Implemented: yes, partially.

Properties that do change:

- aperture;
- velocity after flow recomputation;
- effective connectivity through mesh removal and renumbering.

Properties that do not change explicitly:

- coordinates;
- matrix porosity;
- fracture length.

#### DFN connectivity or topology

Implemented: yes, in an event-based sense.

Mechanism:

- if aperture reaches zero, the mesh is considered closed;
- zero-aperture meshes are removed;
- the network is rebuilt;
- node numbering is rebuilt;
- the backbone is recomputed;
- disconnected particles may be removed from the system.

Where:

- `RebuildNetworkRemovingZeroAperture()`
- `return_backbone()`
- `RebuildNodeNumbering()`
- post-rebuild particle remapping in `Particles_Transport()`

This is not smooth geometric evolution. It is threshold/event-based topological pruning.

### What is only implied

- The intended physical interpretation is reactive dissolution widening fractures, then increasing hydraulic conductance and possibly redirecting transport.
- The code comments suggest geochemical motivation, but the actual implementation is empirical and particle-throughput-based rather than chemistry-balance-based.

### What is missing

- No porosity evolution.
- No mineral mass balance.
- No concentration field.
- No explicit permeability tensor.
- No aperture-dependent geometric widening of coordinates.
- No mechanical deformation.
- No precipitation or clogging model beyond aperture clamping to zero if a formula ever produced it.

## F. Feedback between transport and DFN evolution

The active code implements a genuine feedback loop.

### Forward direction: transport to DFN

At the end of each reaction step:

1. Each active particle contributes advective distance traveled in each mesh.
2. Per mesh, the code estimates an effective passage count:
   - integer full traversals,
   - plus fractional traversals.
3. That effective count `Nb_eff` is normalized by `reference_injected_count_dt`.
4. `ChangeAperture()` updates the aperture of each affected mesh.

### Backward direction: DFN to transport

After aperture updates:

1. `ComputeFlowVelocities()` recomputes the flow field on the modified DFN.
2. `NormalizeMeshDirections()` enforces positive flow orientation by swapping mesh endpoints when needed.
3. If some apertures are zero, the network is rebuilt and pruned.
4. `net_mesh = net_mesh_modified` assigns the updated DFN to the active transport network.
5. The next reaction step therefore uses updated velocities, updated inlet flux weights, and possibly a changed network topology.

### Coupling character

- aperture update: local and periodic;
- flow recomputation: global and periodic;
- topology change: local trigger, global rebuild;
- transport feedback timing: explicit staggered coupling, not fully continuous.

The transport solver uses a fixed DFN during one reaction step and only updates the DFN between reaction steps.

## G. Major assumptions and simplifications

- Fracture geometry is represented as 2D line segments with scalar aperture.
- Flow is solved on the DFN graph only.
- Matrix effects enter transport through diffusion-related residence-time and transfer formulations, not through a resolved matrix concentration field.
- Aperture evolution is empirical and driven by particle traffic, not by concentration or reaction kinetics derived from species balance.
- The same matrix porosity and diffusion coefficient are global constants.
- The DFN is periodically updated at a fixed user-defined `reaction_dt`.
- Flow-induced path selection at intersections is based on a flux-proportional / streamline-routing rule.

## H. Current limitations and missing components

- The coupling variable is particle passage, not solute mass or concentration.
- No explicit reactive chemistry state is solved.
- No porosity update is implemented.
- No separate permeability state is implemented.
- No continuous topology evolution exists; the network changes only when thresholds are crossed.
- No geometric widening of fracture coordinates is implemented.
- The active stepping code uses `Get_Total_Time_From_Advec_Time1D()` in both matrix modes, so the residence-time correction is still the same 1D formula even in the finite-matrix stepping path.
- Some helper routines suggest broader boundary-handling capabilities than are clearly used by the active main path.
- The code still contains older transport routines (`particle_displacement`) and commented alternatives, so the implementation history is visible inside the active source tree.
- There are potential code-quality issues in the active source, for example `fabs(current_mesh.velocity==0.0)` appears to apply `fabs` to a boolean expression rather than to the velocity value itself. This does not change the high-level design, but it matters for robustness.

## I. Text-based workflow diagram

```text
Input/File_names.txt
    ->
Parameters::read_param()
    ->
Domain + selected DFN file + simulation controls
    ->
NetworkMeshes construction
    ->
return_backbone()
    ->
initial DFN flow solution via ComputeFlowVelocities()
    ->
Transport::Particles_Transport()
    ->
Particles_Injection()
    ->
for each reaction step:
    ->
inject ready particles at left-boundary inlet meshes
    ->
for each active particle:
    ->
particle_displacement_step()
    ->
    infinite_matrix_displacement_step()
    or
    finite_matrix_displacement_step()
    ->
    if node reached:
        Mesh_Neighbour_Choice_And_Test()
    ->
    if right boundary reached:
        store arrival time and deactivate particle
    ->
accumulate moved_distances per mesh
    ->
compute Nb_eff per mesh
    ->
ChangeAperture(mesh, reaction_dt, Nb_eff, reference_injected_count_dt)
    ->
if zero aperture or low velocity:
    rebuild DFN and backbone
    remap or remove particles
    ->
ComputeFlowVelocities() on updated DFN
    ->
next reaction step
    ->
Results::post_processing()
    ->
Results::writing()
    ->
WritePositionSnapshotsCSV() + DFN output files
```

## J. Short bullet-point summary

- The active solver is a DFN particle-tracking code with periodic aperture updating.
- Transport is executed on a fixed DFN over each reaction step `reaction_dt`.
- Particle traffic through each fracture mesh is converted into an effective passage measure `Nb_eff`.
- Aperture is updated locally from `Nb_eff` using an empirical law.
- Updated aperture feeds back into hydraulic conductance and transmissivity through the flow recomputation.
- If fractures close or become hydraulically inactive, the DFN is rebuilt and particles may be remapped or removed.
- No dynamic porosity, concentration, or full reactive chemistry model is implemented.

## Important source files for understanding the model

- `Code/src/PERFORM.cpp`
- `Code/src/Transport/Transport.cpp`
- `Code/src/Transport/Transfer.cpp`
- `Code/src/Transport/Projection.cpp`
- `Code/src/Domain_Definition/NetworkMeshes.cpp`
- `Code/src/Domain_Definition/DFNComputation.cpp`
- `Code/src/Domain_Definition/FractureMesh.cpp`
- `Code/src/Chemistry/Chemistry.cpp`
- `Code/src/Input_Output/Parameters.cpp`
- `Code/src/Input_Output/Results.cpp`

## Main functions/subroutines for particle transport

- `Transport::Particles_Transport()`
- `Transport::Particles_Injection()`
- `Transport::particle_displacement_step()`
- `Transport::infinite_matrix_displacement_step()`
- `Transport::finite_matrix_displacement_step()`
- `Transport::Mesh_Neighbour_Choice_And_Test()`
- `Transport::select_neighbour_slr()`
- `Transport::Transfer_Probability_and_Transfer_Time()`
- `Transport::Transfer_Probability_Computation()`
- `Transport::Transfer_Time_Computation()`
- `FractureMesh::Get_Total_Time_From_Advec_Time1D()`
- `FractureMesh::Advection_Time_Computation()`

## Main functions/subroutines for DFN updating

- `NetworkMeshes::ChangeAperture()`
- `Chemistry::ComputeAperture()`
- `RebuildNetworkRemovingZeroAperture()`
- `NetworkMeshes::return_backbone()`
- `RebuildNodeNumbering()`
- `ComputeFlowVelocities()`
- `NetworkMeshes::EvaluateFlowVelocities()`
- `ZeroApertureIfLowVelocity()`

## Short summary of the current coupling strategy

The current coupling strategy is an operator-split, particle-throughput-driven DFN evolution scheme. During each fixed reaction interval, particles move on a frozen DFN. After that interval, each fracture mesh receives a local aperture update based on an effective particle-passage measure accumulated from advective travel distance. The modified apertures are then used to recompute the global DFN flow field, and if some fractures become hydraulically inactive or close, the network is rebuilt before the next transport interval begins.
