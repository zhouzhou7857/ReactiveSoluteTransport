# Chemistry Input Status

The active chemistry path no longer reads a chemistry input file.

`Vref` is now read from the simulation input file as its last optional value.
The older chemistry files in this folder are kept only as legacy records.

## Active Simulation Input Extension

After the existing simulation entries, one extra optional line can be added:

11. `Vref` [m^3]

## Current Active Chemistry Law

The active mineral-volume law is defined in code:

`DeltaV(t) = A1 * (1 - exp(-k1 * t)) + A2 * (1 - exp(-k2 * t)) + L * t`

The simulation input file only sets `Vref`.
It does not set `A1`, `k1`, `A2`, `k2`, or `L`.

## Example

```txt
...
100
1.0e-3
```

This means:

- `Vref = 1.0e-3 m^3`
- the preceding line above it is still the reaction update time step
