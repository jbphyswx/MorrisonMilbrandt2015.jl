# MorrisonMilbrandt2015.jl

Homogeneous-parcel condensation and deposition rates. Three schemes share C2–C7
algebra; they differ in how they treat temperature and events.

| Scheme | Physics |
|--------|---------|
| [`MM2015PiecewiseLinear`](@ref) | `S = δ / (τ Γ)` until the next milestone. External forcing is `A_c` without the WBF addend. |
| [`MM2015FixedT`](@ref) | Appendix C at frozen *T*. Saturation time is the C5 invert. |
| [`MM2015`](@ref) | Appendix C with *T*-updating: frozen-coefficient C6 mass; saturation and freeze are residual roots after a NonEquilibrium `pθq` update. |

## Units and contracts

- Public rates `(S_ql, S_qi)` are specific humidities (kg/kg of moist air per second).
- Freeze events use **`T_triple`**, where `q_sl` and `q_si` cross.
- `Δt == 0` returns `(0, 0)`. `q_tot ≥ 1` is an error.
- A host applies cloud fraction as an area weight of homogeneous sub-regions.
- The C7 liquid-to-ice saturation-frame offset is intrinsic; hosts apply area
  weighting outside the homogeneous kernel.

## Thermodynamics

Kernels take unpacked scalars (`g`, `L`, `c_p`, `∂q*/∂T`, …).
[`MM2015PiecewiseLinear`](@ref) and [`MM2015FixedT`](@ref) use those scalars only.

[`MM2015`](@ref) recovers *T* from `θ_liq_ice` at fixed prognostic `q_liq`, `q_ice`
([`air_temperature_noneq_pθq`](@ref)). Host `dTdt` evolves `θ_li`; with `dTdt = 0`,
`θ_li` is conserved on a segment.

Use [`DefaultThermodynamicsBackend`](@ref), or Thermodynamics.jl methods on
`ThermodynamicsParameters` (params plus variables).

```@contents
```
