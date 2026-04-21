# Chemistry Input Status

The active chemistry path no longer reads a chemistry input file.

`Vref` is now read from the simulation input file, together with fracture thickness.
The older chemistry files in this folder are kept only as legacy records.

## Active Simulation Input Extension

After the existing simulation entries, two active chemistry-related lines can be added:

11. `Vref` [m^3]
12. `fracture thickness` [m]

## Current Active Chemistry Law

The active mineral-volume law is defined in code:

`DeltaV(t) = A1 * (1 - exp(-k1 * t)) + A2 * (1 - exp(-k2 * t)) + L * t`

The simulation input file sets `Vref` and fracture thickness.
It does not set `A1`, `k1`, `A2`, `k2`, or `L`.

## Example

```txt
...
100
1.0e-3
1.0e-1
```

This means:

- `Vref = 1.0e-3 m^3`
- `fracture thickness = 1.0e-1 m`
- the preceding line above them is still the reaction update time step
