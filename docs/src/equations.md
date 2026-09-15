# Appendix C algebra

Morrison & Milbrandt (2015) Appendix C linearizes the mixed-phase supersaturation budget
at a frozen thermodynamic state. Shared implementations live in `src/equations.jl`.

| Symbol | Role |
|--------|------|
| **C2** `τ` | Combined relaxation time ([`τ_func`](@ref) / [`τ_func_combined`](@ref)) |
| **C3** `Γ_l`, `Γ_i` | Psychrometric factors from [`prepare_coefficients`](@ref) |
| **C4** `A_c` | Constant-flux forcing. [`A_c_func_no_WBF`](@ref) is external only; [`A_c_func`](@ref) adds the WBF term |
| **C5** `δ(t)` | `A_c τ + (δ_0 − A_c τ) e^{−t/τ}` ([`δ_func_EP`](@ref)) |
| **C6/C7** | Mean condensation / deposition over a segment ([`S_ql_func`](@ref), [`S_qi_func`](@ref)) |

The C7 liquid-to-ice saturation-frame offset is part of the homogeneous equation.
Subgrid overlap or cloud-fraction weighting belongs outside this kernel.
