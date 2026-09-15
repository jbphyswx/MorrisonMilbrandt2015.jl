# MorrisonMilbrandt2015.jl

Homogeneous-parcel condensation and deposition rates, from Morrison & Milbrandt (2015)
Appendix C (C2–C7) plus a piecewise-linear limiter that uses the same `Γ` and regimes.

| Scheme | Physics |
|--------|---------|
| **MM2015PiecewiseLinear** | `S = δ / (τ Γ)` until the next milestone (condensate out, saturation, or a δ-stationary point). External forcing is `A_c` without the WBF addend. |
| **MM2015FixedT** | Event-driven Appendix C at frozen *T*. Saturation time is the C5 invert; condensate out is Lambert. |
| **MM2015** | Appendix C with *T*-updating: C6 mass at frozen coefficients, then NonEquilibrium `pθq` residuals `F_δ`, `F_δi`, `F_T` for saturation and `T_triple`. |

Public rates `(S_ql, S_qi)` are specific humidities. Freeze events use `T_triple`, where
`q_sl` and `q_si` cross. `Δt == 0` returns `(0, 0)`. `q_tot ≥ 1` is an error.

The C7 liquid-to-ice saturation-frame offset is intrinsic to the homogeneous
equations. A host with cloud fraction applies area weighting outside this kernel.

## Thermodynamics

PiecewiseLinear and EP take unpacked scalars (`g`, `L`, `c_p`, `∂q*/∂T`, …).

MM2015 recovers *T* from `θ_liq_ice` at fixed prognostic `q_liq`, `q_ice`
(NonEquilibrium `pθq`). Host `dTdt` evolves `θ_li`; with `dTdt = 0`, `θ_li` is
conserved on a segment. Use `DefaultThermodynamicsBackend`, or `using Thermodynamics`
for Thermodynamics.jl methods on `ThermodynamicsParameters` (params plus variables).

```julia
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

S_ql, S_qi = MM2015.morrison_milbrandt_2015(MM2015.MM2015FixedT(), inputs)
S_ql, S_qi = MM2015.morrison_milbrandt_2015(MM2015.MM2015(), inputs; thermo)
```
