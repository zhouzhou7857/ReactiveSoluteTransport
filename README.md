# ReactiveSoluteTransport (V5-final)

ReactiveSoluteTransport is a 2D DFN-based reactive solute transport code.
The current V5-final branch is mainly tuned for gypsum precipitation in supersaturated water, and similar cases where reaction volume change is known from a reference chemistry curve.

## What This Version Is Best For

- Supersaturated Ca-SO4 water in gypsum fractures.
- Precipitation-driven aperture reduction and flow redistribution.
- Cases where chemistry can be represented by a reference cumulative mineral-volume curve from PHREEQC.

## Core Chemistry Mechanism (Simple)

This version uses one shared reference reaction curve, then scales it for each particle by the representative injected volume.

1. Build a reference cumulative curve from PHREEQC: `F_ref(t)`.
We used an empirical equation (Rosenberg et al, 2011) that was fitted by a two-stage ational curve : P(t)=P∞​[1−w(1+τ10​t​)−1/9−(1−w)(1+τ2​t​)−1])
In PHREEQC, V_ref is set as the total volume of injection fluid. Surface reaction area is set as the surface area of fracture.

3. For a transport step from `t_start` to `t_end`, compute:

```text
DeltaV_ref = F_ref(t_end) - F_ref(t_start)
```

3. Scale to each particle volume:

```text
DeltaV_particle = DeltaV_ref * (V_particle / V_ref)
```

4. Convert `DeltaV_particle` to aperture change and update DFN.

For gypsum precipitation, `DeltaV` is negative (mineral accumulates in fracture void space), so aperture decreases over time.

## Main Workflow

1. Read domain/simulation/DFN inputs.
2. Build DFN and solve flow. (Pre-processes)
3. Run PHREEQC model, get the cumulative volume function.
4. Run particle transport.
5. Couple chemistry volume change to aperture update.
6. Recompute flow after geometry evolution.
7. Export BTC, particle snapshots, and DFN evolution.

## Latest Benchmark Result Folder

The latest full benchmark case folder in `Output` is:

`Output/benchmark_gypsum_single_fracture_10cm_ap5mm_total_injection_reference_volume_precipitation_20260426_174012/`

This folder includes case index, summaries, validation report, fitted PHREEQC curves, and archived script snapshots.

## Source Code Entry Points

- `Code/src/PERFORM.cpp`
- `Code/src/Transport/Transport.cpp`
- `Code/src/Chemistry/Chemistry.cpp`
- `Code/src/Input_Output/Parameters.cpp`
