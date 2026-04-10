# Decay + Reactivity Cases

These three chemistry cases are designed to satisfy two conditions at the same
time:

1. they all have appreciable decay
2. they have different reactivity strengths

All three use the same transport setup:

- `Domain1.txt`
- `Simulation.txt`
- `DFN_newformat_17.txt`

## Core Design Choice

To keep decay clearly visible, all three cases use the same strong decay rate:

`k_decay = 5.355e-1 1/s`

This is the same order used in the previous `with decay` case, chosen so that
the particle-carried reactant is depleted very quickly relative to the early
crossing time.

To vary reactivity strength without removing the decay effect, only the initial
reactive concentration `C0` is changed.

Since the small-step geometry response scales approximately with:

`effective reactivity ~ C0 * k_decay`

holding `k_decay` fixed and increasing `C0` gives:

- stronger initial geometry impact
- but still rapid reactant exhaustion

## The Three Cases

### `Chemistry_decay_reactivity_case_1.txt`

- `C0 = 0.01`
- `k = 5.355e-1`

This is a mild-with-decay case. It is about:

- `2.83x` stronger than the previous `with decay` baseline (`C0 = 0.00353`)

### `Chemistry_decay_reactivity_case_2.txt`

- `C0 = 0.03`
- `k = 5.355e-1`

This is an intermediate-with-decay case. It is about:

- `8.50x` stronger than the previous `with decay` baseline

### `Chemistry_decay_reactivity_case_3.txt`

- `C0 = 0.1`
- `k = 5.355e-1`

This is a strong-with-decay case. It is about:

- `28.3x` stronger than the previous `with decay` baseline

## Interpretation

These three cases are better suited than the earlier `Chemistry_reactivity_case_*`
files if your goal is specifically:

- to compare different reactivity levels
- while also retaining a clear concentration-depletion effect

The earlier `Chemistry_reactivity_case_*` files increase reactivity mainly in a
quasi-no-decay regime. The new `Chemistry_decay_reactivity_case_*` files are the
correct set for "different reactivity + observable decay" studies.

## Matching File Selectors

- `File_names_decay_reactivity_case_1.txt`
- `File_names_decay_reactivity_case_2.txt`
- `File_names_decay_reactivity_case_3.txt`

## Suggested Run Order

1. previous `with decay` baseline
2. `Chemistry_decay_reactivity_case_1.txt`
3. `Chemistry_decay_reactivity_case_2.txt`
4. `Chemistry_decay_reactivity_case_3.txt`

This run order should make it easier to identify:

- whether increasing reactivity amplifies aperture change
- whether decay still suppresses downstream geometry updates
- when geometry feedback becomes strong enough to alter transport curves
