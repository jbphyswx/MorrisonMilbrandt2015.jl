# MM2015: the parcel model integrated by exprb43, with events located on the step map.

"""Thermodynamic quantities of the parcel model at time `t` and state `u = (x_l, x_i, T)`, per unit basis mass."""
@inline function parcel_state(problem::MM2015Problem{FT}, t::FT, u::Vec3{FT}) where {FT}
    (; basis, thermo, state, forcing) = problem
    x_l, x_i, T = u
    p = state.p + forcing.dpdt * t
    x_t = state.x_tot + forcing.dx_vap_dt * t
    q_d = dry_fraction(basis, x_t)
    q_t = specific_humidity(basis, x_t, q_d)
    q_l = specific_humidity(basis, x_l, q_d)
    q_i = specific_humidity(basis, x_i, q_d)
    liquid = saturation(thermo, T, p, Liquid())
    ice = saturation(thermo, T, p, Ice())
    c_pm = cp_m(thermo, q_t, q_l, q_i)
    R_m = gas_constant_air(thermo, q_t, q_l, q_i)
    c_p = basis_heat_capacity(basis, c_pm, q_d)
    liquid_basis = basis_saturation(basis, liquid, q_d)
    ice_basis = basis_saturation(basis, ice, q_d)
    x_sl, xsl_T, xsl_p, xsl_xt = liquid_basis.x_s, liquid_basis.dx_s_dT, liquid_basis.dx_s_dp, liquid_basis.dx_s_dxt
    x_si, xsi_T, xsi_p, xsi_xt = ice_basis.x_s, ice_basis.dx_s_dT, ice_basis.dx_s_dp, ice_basis.dx_s_dxt
    return (;
        T, p, x_l, x_i, q_d, q_l, q_i, liquid, ice, c_pm, R_m, c_p, ρc_pm = p * c_pm / (R_m * T), x_v = x_t - x_l - x_i,
        x_sl, xsl_T, xsl_p, xsl_xt, x_si, xsi_T, xsi_p, xsi_xt,
        Γ_l = 1 + liquid.L * xsl_T / c_p, Γ_i = 1 + ice.L * xsi_T / c_p,
    )
end

"""Liquid and ice rates and the temperature tendency `(; liq, ice, T)` at the parcel state `s` with active phases `a`."""
@inline function parcel_tendencies(problem::MM2015Problem, a::ActivePhases, s)
    (; τ_liq, τ_ice) = problem.timescales
    (; dpdt, dTdt_external) = problem.forcing
    S_l = a.liquid ? (s.x_v - s.x_sl) / (τ_liq * s.Γ_l) : zero(s.x_v)
    S_i = a.ice ? (s.x_v - s.x_si) / (τ_ice * s.Γ_i) : zero(s.x_v)
    return (; liq = S_l, ice = S_i, T = dTdt_external + dpdt / s.ρc_pm + (s.liquid.L * S_l + s.ice.L * S_i) / s.c_p)
end

"""Right side of the parcel model with fixed active phases, called as `rhs(t, u)`."""
struct ParcelRHS{P}
    problem::P
    active::ActivePhases
end

@inline (rhs::ParcelRHS)(t, u) = parcel_tendencies(rhs.problem, rhs.active, parcel_state(rhs.problem, t, u))

"""
    parcel_linearization(problem, a, s) -> (; F, J, v)

Right side `F`, Jacobian `J = ∂F/∂u` (by columns), and `v = ∂F/∂t` of the parcel model with active
phases `a` at the [`parcel_state`](@ref) `s`.
"""
@inline function parcel_linearization(problem::MM2015Problem{FT}, a::ActivePhases, s) where {FT}
    (; basis, thermo, timescales, forcing) = problem
    F = parcel_tendencies(problem, a, s)
    S_l, S_i = F.liq, F.ice
    ṗ, F_t = forcing.dpdt, forcing.dx_vap_dt
    c_v, c_l, c_i, c_d = cp_v(thermo, FT), cp_l(thermo, FT), cp_i(thermo, FT), cp_d(thermo, FT)
    Rv, Rd = R_v(thermo, FT), R_d(thermo, FT)
    τ_l, τ_i = timescales.τ_liq, timescales.τ_ice
    L_v, L_s = s.liquid.L, s.ice.L
    c_p, Q = s.c_p, s.ρc_pm

    liquid_second = basis_saturation_second(basis, s.liquid, s.q_d)
    ice_second = basis_saturation_second(basis, s.ice, s.q_d)
    xsl_TT, xsl_Tp, xsl_Txt = liquid_second.x_s_TT, liquid_second.x_s_Tp, liquid_second.x_s_Txt
    xsi_TT, xsi_Tp, xsi_Txt = ice_second.x_s_TT, ice_second.x_s_Tp, ice_second.x_s_Txt
    (; dq_t, dq_l, dq_i) = total_water_derivatives(basis, s.q_d, s.q_l, s.q_i)
    dc_pm_dxt = (c_v - c_d) * dq_t + (c_l - c_v) * dq_l + (c_i - c_v) * dq_i
    dR_m_dxt = (Rv - Rd) * dq_t - Rv * (dq_l + dq_i)
    dc_p_dxt = basis_heat_capacity_derivative(basis, dc_pm_dxt, s.c_pm, s.q_d)
    dc_p_dxl, dc_p_dxi = c_l - c_v, c_i - c_v
    w_q = specific_humidity(basis, one(FT), s.q_d)
    dQ_dxl = Q * w_q * (dc_p_dxl / s.c_pm + Rv / s.R_m)
    dQ_dxi = Q * w_q * (dc_p_dxi / s.c_pm + Rv / s.R_m)
    dQ_dT = -Q / s.T
    dQ_dt = Q * (ṗ / s.p + F_t * (dc_pm_dxt / s.c_pm - dR_m_dxt / s.R_m))

    # rate S = δ / (τ Γ): ∂S = (∂δ − S τ ∂Γ) / (τ Γ)
    dSl = if a.liquid
        Γ = s.Γ_l
        g = τ_l * Γ
        dΓ_xl, dΓ_xi = -(Γ - 1) * dc_p_dxl / c_p, -(Γ - 1) * dc_p_dxi / c_p
        dΓ_T = (s.liquid.dL_dT * s.xsl_T + L_v * xsl_TT) / c_p
        dΓ_t = (L_v * (xsl_Tp * ṗ + xsl_Txt * F_t) - (Γ - 1) * dc_p_dxt * F_t) / c_p
        dδ_t = F_t * (1 - s.xsl_xt) - s.xsl_p * ṗ
        ((-1 - S_l * τ_l * dΓ_xl) / g, (-1 - S_l * τ_l * dΓ_xi) / g, (-s.xsl_T - S_l * τ_l * dΓ_T) / g, (dδ_t - S_l * τ_l * dΓ_t) / g)
    else
        (zero(FT), zero(FT), zero(FT), zero(FT))
    end
    dSi = if a.ice
        Γ = s.Γ_i
        g = τ_i * Γ
        dΓ_xl, dΓ_xi = -(Γ - 1) * dc_p_dxl / c_p, -(Γ - 1) * dc_p_dxi / c_p
        dΓ_T = (s.ice.dL_dT * s.xsi_T + L_s * xsi_TT) / c_p
        dΓ_t = (L_s * (xsi_Tp * ṗ + xsi_Txt * F_t) - (Γ - 1) * dc_p_dxt * F_t) / c_p
        dδ_t = F_t * (1 - s.xsi_xt) - s.xsi_p * ṗ
        ((-1 - S_i * τ_i * dΓ_xl) / g, (-1 - S_i * τ_i * dΓ_xi) / g, (-s.xsi_T - S_i * τ_i * dΓ_T) / g, (dδ_t - S_i * τ_i * dΓ_t) / g)
    else
        (zero(FT), zero(FT), zero(FT), zero(FT))
    end
    H = L_v * S_l + L_s * S_i
    dT_xl = -ṗ * dQ_dxl / Q^2 + (L_v * dSl[1] + L_s * dSi[1]) / c_p - H * dc_p_dxl / c_p^2
    dT_xi = -ṗ * dQ_dxi / Q^2 + (L_v * dSl[2] + L_s * dSi[2]) / c_p - H * dc_p_dxi / c_p^2
    dT_T = -ṗ * dQ_dT / Q^2 + (s.liquid.dL_dT * S_l + L_v * dSl[3] + s.ice.dL_dT * S_i + L_s * dSi[3]) / c_p
    dT_t = -ṗ * dQ_dt / Q^2 + (L_v * dSl[4] + L_s * dSi[4]) / c_p - H * dc_p_dxt * F_t / c_p^2
    J = (dSl[1], dSi[1], dT_xl, dSl[2], dSi[2], dT_xi, dSl[3], dSi[3], dT_T)
    v = (dSl[4], dSi[4], dT_t)
    return (; F, J, v)
end

"""Time derivatives `(; dδ, dδ_i)` of `δ` and `δ_i` at the [`parcel_state`](@ref) `s`, from [`parcel_tendencies`](@ref)."""
@inline function supersaturation_rates(problem::MM2015Problem, s, tendencies)
    S_l, S_i, dT = tendencies.liq, tendencies.ice, tendencies.T
    (; dpdt, dx_vap_dt) = problem.forcing
    dx_v = dx_vap_dt - S_l - S_i
    return (;
        dδ = dx_v - s.xsl_T * dT - s.xsl_p * dpdt - s.xsl_xt * dx_vap_dt,
        dδ_i = dx_v - s.xsi_T * dT - s.xsi_p * dpdt - s.xsi_xt * dx_vap_dt,
    )
end

"""Whether the parcel at the [`parcel_state`](@ref) `s` counts as below the triple point."""
@inline function below_triple(problem::MM2015Problem{FT}, s) where {FT}
    T_tr = T_triple(problem.thermo, FT)
    s.T == T_tr || return s.T < T_tr
    (; dpdt, dTdt_external) = problem.forcing
    return dTdt_external + dpdt / s.ρc_pm < 0
end

"""
    active_phases(problem, s, below, thresholds = Thresholds(0, 0, 0)) -> ActivePhases

The rule of [`active_phases`](@ref active_phases(::Coefficients, ::Any, ::Any, ::Any, ::Any)) at the
[`parcel_state`](@ref) `s`, with the time derivatives of the parcel model on the saturation boundaries.
A phase without mass forms only when its supersaturation exceeds `thresholds.δ_min` (liquid) or
`thresholds.δ_i_min` (ice), or lies within it, is not negative, and rises.
"""
@inline function active_phases(
    problem::MM2015Problem{FT},
    s,
    below::Bool,
    thresholds::Thresholds{FT} = Thresholds{FT}(zero(FT), zero(FT), zero(FT)),
) where {FT}
    (; x_l, x_i) = s
    δ, δ_i = s.x_v - s.x_sl, s.x_v - s.x_si
    rounded = rounded_supersaturations(s, thresholds)
    δ_form, δ_i_form = rounded.δ, rounded.δ_i
    ice_interior = below ? (x_i > 0 || δ_i_form > 0) : (x_i > 0 && δ_i < 0)
    liquid = x_l > 0 || δ_form > 0 ||
             (iszero(x_l) && iszero(δ_form) && δ ≥ 0 &&
              supersaturation_rates(problem, s, parcel_tendencies(problem, ActivePhases(false, ice_interior), s)).dδ > 0)
    ice_saturation_rate() = supersaturation_rates(problem, s, parcel_tendencies(problem, ActivePhases(liquid, false), s)).dδ_i
    ice = if below
        ice_interior || (iszero(x_i) && iszero(δ_i_form) && δ_i ≥ 0 && ice_saturation_rate() > 0)
    else
        ice_interior || (x_i > 0 && iszero(δ_i) && ice_saturation_rate() < 0)
    end
    return ActivePhases(liquid, ice)
end

"""Supersaturations `(; δ, δ_i)` over liquid and over ice at the [`parcel_state`](@ref) `s`, zero below `thresholds.δ_min` and `thresholds.δ_i_min` in magnitude."""
@inline function rounded_supersaturations(s, thresholds::Thresholds)
    δ, δ_i = s.x_v - s.x_sl, s.x_v - s.x_si
    return (; δ = abs(δ) < thresholds.δ_min ? zero(δ) : δ, δ_i = abs(δ_i) < thresholds.δ_i_min ? zero(δ_i) : δ_i)
end

"""
    event_indicators(problem, a, below, s, thresholds = Thresholds(0, 0, 0)) -> NTuple{4}

Values at the [`parcel_state`](@ref) `s` that end a segment when one becomes positive: `-x_l`
(active liquid) or `δ`; `-x_i` for active ice; `δ_i` for inactive ice below the triple point and
active ice above it, `-δ_i` for inactive ice with mass above it; and the signed distance past
`T_triple`. `-Inf` marks an unused indicator. For a phase without mass the supersaturation is that of
[`rounded_supersaturations`](@ref).
"""
@inline function event_indicators(
    problem::MM2015Problem{FT},
    a::ActivePhases,
    below::Bool,
    s,
    thresholds::Thresholds{FT} = Thresholds{FT}(zero(FT), zero(FT), zero(FT)),
) where {FT}
    none = -FT(Inf)
    (; x_l, x_i, T) = s
    δ_i = s.x_v - s.x_si
    rounded = rounded_supersaturations(s, thresholds)
    δ_form, δ_i_form = rounded.δ, rounded.δ_i
    T_tr = T_triple(problem.thermo, FT)
    liquid = a.liquid ? -x_l : δ_form
    ice_mass = a.ice ? -x_i : none
    ice_saturation = below ? (a.ice ? none : δ_i_form) : (a.ice ? δ_i : (x_i > 0 ? -δ_i : none))
    return (liquid, ice_mass, ice_saturation, below ? T - T_tr : T_tr - T)
end

"""Time derivatives of [`event_indicators`](@ref) along the parcel model at the [`parcel_state`](@ref) `s`."""
@inline function event_indicator_rates(problem::MM2015Problem{FT}, a::ActivePhases, below::Bool, s) where {FT}
    tendencies = parcel_tendencies(problem, a, s)
    S_l, S_i, dT = tendencies.liq, tendencies.ice, tendencies.T
    (; dδ, dδ_i) = supersaturation_rates(problem, s, tendencies)
    x_i = s.x_i
    liquid = a.liquid ? -S_l : dδ
    ice_mass = a.ice ? -S_i : zero(FT)
    ice_saturation = below ? (a.ice ? zero(FT) : dδ_i) : (a.ice ? dδ_i : (x_i > 0 ? -dδ_i : zero(FT)))
    return (liquid, ice_mass, ice_saturation, below ? dT : -dT)
end

"""Largest value on `[0, 1]` of the cubic Hermite interpolant of `(g₀, h ġ₀)` at 0 and `(g₁, h ġ₁)` at 1."""
@inline function hermite_maximum(g₀::FT, g₁::FT, d₀::FT, d₁::FT) where {FT}
    # p(σ) = g₀ + d₀σ + (3(g₁ − g₀) − 2d₀ − d₁)σ² + (2(g₀ − g₁) + d₀ + d₁)σ³
    b, c = 3 * (g₁ - g₀) - 2 * d₀ - d₁, 2 * (g₀ - g₁) + d₀ + d₁
    p(σ) = g₀ + σ * (d₀ + σ * (b + σ * c))
    best = max(g₀, g₁)
    # p′(σ) = d₀ + 2bσ + 3cσ²
    if iszero(c)
        iszero(b) || (σ = -d₀ / (2b); zero(FT) < σ < one(FT) && (best = max(best, p(σ))))
    else
        disc = b^2 - 3 * c * d₀
        if disc ≥ 0
            r = sqrt(disc)
            for σ in ((-b - r) / (3c), (-b + r) / (3c))
                zero(FT) < σ < one(FT) && (best = max(best, p(σ)))
            end
        end
    end
    return best
end

"""`(s, e)` with `s = fl(a + b)` and `a + b = s + e` exactly (Knuth's TwoSum)."""
@inline function two_sum(a::FT, b::FT) where {FT}
    s = a + b
    b′ = s - a
    return s, (a - (s - b′)) + (b - b′)
end

"""The state `u + u_lo` advanced by `Δu`, as a rounded state and its rounding error."""
@inline function compensated_add(u::Vec3{FT}, u_lo::Vec3{FT}, Δu::Vec3{FT}) where {FT}
    sums = map(two_sum, u, Δu .+ u_lo)
    return map(first, sums), map(last, sums)
end

"""The indicators `g` with those not `selected` replaced by `-Inf`."""
@inline masked(g::NTuple{4, FT}, selected::NTuple{4, Bool}) where {FT} = ntuple(k -> selected[k] ? g[k] : -FT(Inf), Val(4))

"""Selection of every indicator except `k`."""
@inline all_but(k::Int) = ntuple(m -> m != k, Val(4))

"""
Event indicators after a step of length `s` from `(t, u + u_lo)` with the linearization `(F, J, v)`,
called as `indicators(s)`; `indicators(s, k)` is the `k`-th.
"""
struct StepIndicators{P, FT}
    problem::P
    active::ActivePhases
    below::Bool
    thresholds::Thresholds{FT}
    t::FT
    u::Vec3{FT}
    u_lo::Vec3{FT}
    F::Vec3{FT}
    J::Mat3{FT}
    v::Vec3{FT}
    scaling::NTuple{2, Mat3{FT}}
end

@inline function (m::StepIndicators)(s)
    Δu, _ = exprb43_step(ParcelRHS(m.problem, m.active), m.t, m.u, s, Tuple(m.F), m.J, m.v, m.scaling)
    state = parcel_state(m.problem, m.t + s, first(compensated_add(m.u, m.u_lo, Δu)))
    return event_indicators(m.problem, m.active, m.below, state, m.thresholds)
end

@inline (m::StepIndicators)(s, k::Int) = m(s)[k]

"""Indicator `k` of `indicators::StepIndicators` as a function of the step length."""
struct StepIndicator{M}
    indicators::M
    k::Int
end

@inline (m::StepIndicator)(s) = m.indicators(s, m.k)

"""
Rounding error of event indicator `k` of [`event_indicators`](@ref) with active phases `a` at the
[`parcel_state`](@ref) `s`: `eps` times the size of the quantities it is formed from.
"""
@inline function indicator_rounding(k::Int, a::ActivePhases, s)
    ε = eps(typeof(s.T))
    k == 1 && return ε * (a.liquid ? abs(s.x_l) : s.x_v + abs(s.xsl_T) * s.T)
    k == 2 && return ε * abs(s.x_i)
    k == 3 && return ε * (s.x_v + abs(s.xsi_T) * s.T)
    return ε * s.T
end

"""Index of the fired indicator (`g₁ > 0`) whose linear interpolation from `g₀` crosses zero first."""
@inline function earliest_fired(g₀::NTuple{4, FT}, g₁::NTuple{4, FT}) where {FT}
    first_k, first_σ = 0, FT(Inf)
    for k in 1:4
        g₁[k] > 0 || continue
        σ = g₀[k] / (g₀[k] - g₁[k])
        σ < first_σ && ((first_k, first_σ) = (k, σ))
    end
    return first_k
end

"""Largest component of the local error estimate relative to `atol + rtol max(|u|, |u_next|)`."""
@inline function scaled_error(estimate::Vec3, u::Vec3, u_next::Vec3, atol::Vec3, rtol)
    return maximum(ntuple(c -> abs(estimate[c]) / (atol[c] + rtol * max(abs(u[c]), abs(u_next[c]))), Val(3)))
end

"""
Whether an indicator can cross zero inside a step of length `h` without a sign change at its ends:
the cubic Hermite interpolant of its values `g`, `g_next` and rates `ġ`, `ġ_next` rises above zero.
"""
@inline function hidden_crossing(g::NTuple{4}, g_next::NTuple{4}, ġ::NTuple{4}, ġ_next::NTuple{4}, h)
    return any(ntuple(k -> isfinite(g[k]) && hermite_maximum(g[k], g_next[k], h * ġ[k], h * ġ_next[k]) > 0, Val(4)))
end

"""Event kind of indicator `k` of [`event_indicators`](@ref) with active phases `a`."""
@inline event_kind(k::Int, a::ActivePhases) =
    k == 1 ? (a.liquid ? LiquidExhausted : LiquidSaturation) : k == 2 ? IceExhausted : k == 3 ? IceSaturation : TripleCrossing

"""
    return_small_condensate(problem, t, u, u_lo, x_min) -> (; u, u_lo, r_l, r_i)

The state `u + u_lo = (x_l, x_i, T)` with liquid and ice below `x_min` returned to the vapor, the
temperature lowered by their latent heat, and the returned amounts `r_l` and `r_i`.
"""
@inline function return_small_condensate(problem::MM2015Problem{FT}, t::FT, u::Vec3{FT}, u_lo::Vec3{FT}, x_min::FT) where {FT}
    x_l, x_i, T = u
    small_l, small_i = zero(FT) < x_l < x_min, zero(FT) < x_i < x_min
    small_l || small_i || return (; u, u_lo, r_l = zero(FT), r_i = zero(FT))
    r_l = small_l ? x_l + u_lo[1] : zero(FT)
    r_i = small_i ? x_i + u_lo[2] : zero(FT)
    s = parcel_state(problem, t, u)
    u = (small_l ? zero(FT) : x_l, small_i ? zero(FT) : x_i, T - (s.liquid.L * r_l + s.ice.L * r_i) / s.c_p)
    u_lo = (small_l ? zero(FT) : u_lo[1], small_i ? zero(FT) : u_lo[2], u_lo[3])
    return (; u, u_lo, r_l, r_i)
end

"""
    step_start(scheme, problem, t, u, u_lo, scaling) -> NamedTuple

Start of an integrator step from `(t, u + u_lo)`: the state with the condensate below the threshold of
`scheme` returned to the vapor, the amounts returned, the parcel state, the active phases, the event
indicators and their rates, and the linearization with its balancing warm-started from `scaling`.
Returns the fields `u`, `u_lo`, `r_l`, `r_i`, `state`, `below`, `a`, `g`, `ġ`, `F`, `J`, `v`, and `scaling`.
"""
@inline function step_start(scheme::MM2015, problem::MM2015Problem{FT}, t::FT, u::Vec3{FT}, u_lo::Vec3{FT}, scaling) where {FT}
    (; thresholds) = scheme
    (; u, u_lo, r_l, r_i) = return_small_condensate(problem, t, u, u_lo, thresholds.x_min)
    state = parcel_state(problem, t, u)
    below = below_triple(problem, state)
    a = active_phases(problem, state, below, thresholds)
    (; F, J, v) = parcel_linearization(problem, a, state)
    g = event_indicators(problem, a, below, state, thresholds)
    ġ = event_indicator_rates(problem, a, below, state)
    scaling = balancing(J, scaling)
    return (; u, u_lo, r_l, r_i, state, below, a, g, ġ, F, J, v, scaling)
end

"""
    evolve(scheme::MM2015, problem, Δt, recorder) -> (; Δx_l, Δx_i)

Liquid and ice increments over `[0, Δt]` of the parcel model, integrated by exprb43 with step-size
control to the tolerances of `scheme`. When event indicators turn positive in an accepted step, the
first to cross is located on the step map by `scheme.root_finder`, and the step is shortened to it.
Each step is passed to `recorder`.
"""
function evolve(scheme::MM2015, problem::MM2015Problem{FT}, Δt::FT, recorder) where {FT}
    rtol, atol_q, atol_T = FT(scheme.rtol), FT(scheme.atol_q), FT(scheme.atol_T)
    rtol ≥ 4 * eps(FT) || throw_rtol_below_precision(scheme.rtol, FT)
    atol = (atol_q, atol_q, atol_T)
    (; thresholds) = scheme
    t = zero(FT)
    (; u, u_lo, r_l, r_i, state, below, a, g, ġ, F, J, v, scaling) = step_start(
        scheme, problem, t, (problem.state.x_liq, problem.state.x_ice, problem.state.T), (zero(FT), zero(FT), zero(FT)),
        identity_scaling(FT))
    Δx_l, Δx_i = -r_l, -r_i
    h = Δt
    for _ in 1:scheme.max_steps
        h = min(h, Δt - t)
        rhs = ParcelRHS(problem, a)
        Δu, error_estimate = exprb43_step(rhs, t, u, h, Tuple(F), J, v, scaling)
        u_next, lo_next = compensated_add(u, u_lo, Δu)
        err = scaled_error(error_estimate, u, u_next, atol, rtol)
        if !(err ≤ 1)
            h *= isfinite(err) ? max(FT(0.2), FT(0.9) / sqrt(sqrt(err))) : FT(0.2)
            continue
        end
        δ, δ_i = state.x_v - state.x_sl, state.x_v - state.x_si
        state_next = parcel_state(problem, t + h, u_next)
        g_next = event_indicators(problem, a, below, state_next, thresholds)
        if any(>(0), g_next)
            indicators = StepIndicators(problem, a, below, thresholds, t, u, u_lo, Tuple(F), J, v, scaling)
            k, s, g_event = earliest_fired(g, g_next), h, g_next
            while true
                # no finer in time than the indicator resolves at its mean rate over (0, s]
                tolerance = max(16 * eps(FT) * (t + h), 8 * indicator_rounding(k, a, state_next) * s / (g_event[k] - g[k]))
                s_lo, s = bracket_root(scheme.root_finder, StepIndicator(indicators, k), zero(FT), s, g[k], g_event[k], tolerance)
                g_event = indicators(s)
                j = earliest_fired(g, masked(g_event, all_but(k)))
                iszero(j) && break
                g_lo = indicators(s_lo)
                # j crossed inside (s_lo, s] with k: both events happen at s
                g_lo[j] > 0 || break
                k, s, g_event = j, s_lo, g_lo
            end
            Δu, _ = exprb43_step(rhs, t, u, s, Tuple(F), J, v, scaling)
            u_event, lo_event = compensated_add(u, u_lo, Δu)
            ġ_event = event_indicator_rates(problem, a, below, parcel_state(problem, t + s, u_event))
            unfired = map(≤(0), g_event)
            if hidden_crossing(masked(g, unfired), masked(g_event, unfired), ġ, ġ_event, s)
                h = s / 2
                continue
            end
            kind = EndOfStep
            liquid_exhausted = ice_exhausted = false
            for m in 1:4
                g_event[m] > 0 || continue
                kind = event_kind(m, a)
                liquid_exhausted |= kind == LiquidExhausted
                ice_exhausted |= kind == IceExhausted
            end
            record!(recorder, t, s, a, u[1], u[2], u[3], δ, δ_i, kind)
            Δx_l += liquid_exhausted ? -(u[1] + u_lo[1]) : Δu[1]
            Δx_i += ice_exhausted ? -(u[2] + u_lo[2]) : Δu[2]
            t += s
            t < Δt || return (; Δx_l, Δx_i)
            u_new = (liquid_exhausted ? zero(FT) : u_event[1], ice_exhausted ? zero(FT) : u_event[2], u_event[3])
            lo_new = (liquid_exhausted ? zero(FT) : lo_event[1], ice_exhausted ? zero(FT) : lo_event[2], lo_event[3])
            (; u, u_lo, r_l, r_i, state, below, a, g, ġ, F, J, v, scaling) = step_start(scheme, problem, t, u_new, lo_new, scaling)
            Δx_l -= r_l
            Δx_i -= r_i
        else
            ġ_next = event_indicator_rates(problem, a, below, state_next)
            if hidden_crossing(g, g_next, ġ, ġ_next, h)
                h /= 2
                continue
            end
            record!(recorder, t, h, a, u[1], u[2], u[3], δ, δ_i, EndOfStep)
            Δx_l += Δu[1]
            Δx_i += Δu[2]
            t += h
            t < Δt || return (; Δx_l, Δx_i)
            if zero(FT) < u_next[1] < thresholds.x_min || zero(FT) < u_next[2] < thresholds.x_min
                (; u, u_lo, r_l, r_i, state, below, a, g, ġ, F, J, v, scaling) = step_start(scheme, problem, t, u_next, lo_next, scaling)
                Δx_l -= r_l
                Δx_i -= r_i
            else
                u, u_lo = u_next, lo_next
                state, g, ġ = state_next, g_next, ġ_next
                (; F, J, v) = parcel_linearization(problem, a, state)
                scaling = balancing(J, scaling)
            end
        end
        h *= min(FT(5), max(FT(0.2), FT(0.9) / sqrt(sqrt(max(err, eps(FT))))))
    end
    throw_max_steps(scheme.max_steps, Δt - t)
end

@noinline throw_rtol_below_precision(rtol, FT) =
    throw(ArgumentError("rtol = $rtol is below the precision of a $FT problem; use MM2015{$FT}() or a larger rtol"))

@noinline throw_max_steps(max_steps, remaining) =
    throw(ErrorException("MorrisonMilbrandt2015: MM2015 exceeded max_steps = $max_steps with $remaining s left"))

"""
    state_at(trajectory::Trajectory{FT, <:MM2015}, t) -> (; δ, δ_i, x_l, x_i, T)

Supersaturation over liquid and over ice, liquid, ice, and temperature at time `t` in `[0, Δt]`, on
the step map of the integrator step that contains `t`.
"""
function state_at(trajectory::Trajectory{FT, <:MM2015}, t::Real) where {FT}
    (; context, segments) = trajectory
    segment = segments[something(findlast(r -> r.t ≤ t, segments), firstindex(segments))]
    u = (segment.x_l, segment.x_i, segment.T)
    (; F, J, v) = parcel_linearization(context, segment.active, parcel_state(context, segment.t, u))
    Δu, _ = exprb43_step(ParcelRHS(context, segment.active), segment.t, u, FT(t) - segment.t, Tuple(F), J, v, balancing(J))
    x_l, x_i, T = u .+ Δu
    s = parcel_state(context, FT(t), (x_l, x_i, T))
    return (; δ = s.x_v - s.x_sl, δ_i = s.x_v - s.x_si, x_l, x_i, T)
end
