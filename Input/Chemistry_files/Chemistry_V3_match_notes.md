# V3 Matching Chemistry Cases

These chemistry files target the current concentration-based DFN-PT-V3 code path
for the case:

- `Domain1.txt`
- `Simulation.txt`
- `DFN_newformat_17.txt`

File format:

1. initial reactive concentration `C0` [mol/m^3]
2. first-order concentration decay `k_decay` [1/s]
3. reactive-to-mineral stoichiometric factor
4. mineral molar volume [m^3/mol]
5. fracture out-of-plane thickness [m]

Calibration logic:

- The old V3 particle-count coupling used an empirical aperture-change rate of
  `1.5e-10 m/s` when the normalized particle-throughput ratio was one.
- For the present concentration-based implementation, the effective small-step
  aperture rate scales as `C0 * k_decay * V_particle * stoich * V_m / L`.
- Using the old V3 output (`Output/0324/V3`) for this DFN, the first non-zero
  arrival time is about `12.8989 s`; the decay case therefore targets practical
  depletion by `2/3` of that time, about `8.60 s`.

Case meaning:

- `Chemistry_V3_match_no_decay.txt`
  keeps `k_decay` very small so concentration loss over one network crossing is
  negligible, while `C0 * k_decay` is scaled to remain close to the old V3
  aperture-growth strength.

- `Chemistry_V3_match_decay.txt`
  uses the same initial aperture-growth strength but with strong depletion,
  parameterized so the carried reactant drops to about 1 percent of its inlet
  value by roughly `8.60 s`.
