# Stronger Reactivity Cases

These three chemistry cases are intended for direct comparison of stronger
particle-driven reactivity on fracture-geometry evolution while keeping the
transport setup fixed:

- `Domain1.txt`
- `Simulation.txt`
- `DFN_newformat_17.txt`

The only change between the three cases is the chemistry file.

## File List

- `Chemistry_reactivity_case_1.txt`
- `Chemistry_reactivity_case_2.txt`
- `Chemistry_reactivity_case_3.txt`

Matching `File_names` helpers:

- `File_names_reactivity_case_1.txt`
- `File_names_reactivity_case_2.txt`
- `File_names_reactivity_case_3.txt`

## Design Logic

In the current code, the small-step aperture response scales approximately with:

`effective reactivity ~ C0 * k_decay`

where:

- `C0` is the initial reactive concentration carried by each particle
- `k_decay` is the first-order concentration decay constant

To make the comparison cleaner, these three cases keep:

- the same DFN
- the same domain and hydraulic forcing
- the same simulation parameters
- the same `k_decay = 1.0e-5`
- the same stoichiometry and mineral-volume conversion

Only `C0` is increased, so the chemistry becomes stronger without changing the
form of the decay law.

## Relative Strength

Relative to `Chemistry_V3_match_no_decay.txt` (`C0 = 189`, `k = 1e-5`):

- `case_1`: `C0 = 500`
  - about `2.65x` the baseline effective reactivity
- `case_2`: `C0 = 2000`
  - about `10.6x` the baseline effective reactivity
- `case_3`: `C0 = 5000`
  - about `26.5x` the baseline effective reactivity

## Suggested Usage

If you want a monotonic geometry comparison, run:

1. baseline no-decay-like case: `Chemistry_V3_match_no_decay.txt`
2. `Chemistry_reactivity_case_1.txt`
3. `Chemistry_reactivity_case_2.txt`
4. `Chemistry_reactivity_case_3.txt`

This sequence should make it easier to identify how increasing reactivity
changes:

- the magnitude of aperture change
- the spatial distribution of aperture change
- the extent to which geometry feedback starts to alter transport
