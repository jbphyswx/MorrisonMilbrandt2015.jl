# Test-only reference solutions of the parcel model: Radau IIA of order 5, events located on the step map.

module ParcelReference

using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

"""Active phases of a segment."""
struct Active
    liquid::Bool
    ice::Bool
end

"""
    FrozenLinear(problem)

Frozen-coefficient dynamics of `MM2015FixedT`: the state `(δ, x_l, x_i)` under the coefficients of
`problem` at its initial state.
"""
struct FrozenLinear{FT}
    k::MM2015.Coefficients{FT}
end
FrozenLinear(problem::MM2015.MM2015Problem) = FrozenLinear(MM2015.coefficients(problem))

"""
    Nonlinear(problem)

The parcel model in the moisture basis of `problem`: the state `(x_l, x_i, T)`, with pressure and
total water linear in time.
"""
struct Nonlinear{P <: MM2015.MM2015Problem}
    problem::P
end

initial_state(m::FrozenLinear) = [m.k.δ, m.k.x_l, m.k.x_i]
initial_state(m::Nonlinear) = (s = m.problem.state; [s.x_liq, s.x_ice, s.T])

"""Liquid and ice of the state `u`."""
condensate(::FrozenLinear, u) = (u[2], u[3])
condensate(::Nonlinear, u) = (u[1], u[2])

condensate_index(::FrozenLinear, liquid::Bool) = liquid ? 2 : 3
condensate_index(::Nonlinear, liquid::Bool) = liquid ? 1 : 2

"""Absolute tolerance and finite-difference scale of each state component."""
component_scales(::FrozenLinear, atol_q, _) = ((atol_q, atol_q, atol_q), (1e-6, 1e-6, 1e-6))
component_scales(::Nonlinear, atol_q, atol_T) = ((atol_q, atol_q, atol_T), (1e-6, 1e-6, 1.0))

supersaturations(m::FrozenLinear, _, u) = (u[1], u[1] + m.k.Δ)
below_triple(m::FrozenLinear, _, _) = m.k.below_triple
triple_value(::FrozenLinear, _, _, u) = -eltype(u)(Inf)

function rhs(m::FrozenLinear, a::Active, _, u)
    k = m.k
    δ = u[1]
    S_l = a.liquid ? δ / (k.τ_l * k.Γ_l) : zero(δ)
    S_i = a.ice ? (δ + k.Δ) / (k.τ_i * k.Γ_i) : zero(δ)
    return [k.A_l - k.Γ_l * S_l - k.α * S_i, S_l, S_i]
end

"""Thermodynamic quantities of the parcel model at time `t` and state `u`, per unit basis mass."""
function parcel(m::Nonlinear, t, u)
    (; basis, thermo, state, forcing) = m.problem
    p = state.p + forcing.dpdt * t
    x_t = state.x_tot + forcing.dx_vap_dt * t
    return parcel(basis, thermo, u[3], p, x_t, u[1], u[2])
end

function parcel(::MM2015.SpecificHumidity, thermo, T, p, q_t, q_l, q_i)
    q_d = 1 - q_t
    liq = MM2015.saturation(thermo, T, p, MM2015.Liquid())
    ice = MM2015.saturation(thermo, T, p, MM2015.Ice())
    c_pm = MM2015.cp_m(thermo, q_t, q_l, q_i)
    ρ = MM2015.air_density(thermo, T, p, q_t, q_l, q_i)
    return (;
        x_v = q_t - q_l - q_i, x_sl = q_d * liq.r, x_si = q_d * ice.r,
        Γ_l = 1 + liq.L / c_pm * q_d * liq.dr_dT, Γ_i = 1 + ice.L / c_pm * q_d * ice.dr_dT,
        L_v = liq.L, L_s = ice.L, c_p = c_pm, ρc_p = ρ * c_pm,
    )
end

function parcel(::MM2015.DryAirMixingRatio, thermo, T, p, r_t, r_l, r_i)
    q_d = 1 / (1 + r_t)
    liq = MM2015.saturation(thermo, T, p, MM2015.Liquid())
    ice = MM2015.saturation(thermo, T, p, MM2015.Ice())
    c_pm = MM2015.cp_m(thermo, q_d * r_t, q_d * r_l, q_d * r_i)
    ρ = MM2015.air_density(thermo, T, p, q_d * r_t, q_d * r_l, q_d * r_i)
    c_p = c_pm / q_d
    return (;
        x_v = r_t - r_l - r_i, x_sl = liq.r, x_si = ice.r,
        Γ_l = 1 + liq.L / c_p * liq.dr_dT, Γ_i = 1 + ice.L / c_p * ice.dr_dT,
        L_v = liq.L, L_s = ice.L, c_p, ρc_p = ρ * c_pm,
    )
end

function supersaturations(m::Nonlinear, t, u)
    s = parcel(m, t, u)
    return (s.x_v - s.x_sl, s.x_v - s.x_si)
end

function rhs(m::Nonlinear, a::Active, t, u)
    s = parcel(m, t, u)
    (; τ_liq, τ_ice) = m.problem.timescales
    (; dpdt, dTdt_external) = m.problem.forcing
    S_l = a.liquid ? (s.x_v - s.x_sl) / (τ_liq * s.Γ_l) : zero(s.x_v)
    S_i = a.ice ? (s.x_v - s.x_si) / (τ_ice * s.Γ_i) : zero(s.x_v)
    return [S_l, S_i, dTdt_external + dpdt / s.ρc_p + (s.L_v * S_l + s.L_s * S_i) / s.c_p]
end

"""Temperature tendency without latent heating [K s⁻¹]."""
function non_latent_tendency(m::Nonlinear, t, u)
    (; dpdt, dTdt_external) = m.problem.forcing
    return dTdt_external + dpdt / parcel(m, t, u).ρc_p
end

function below_triple(m::Nonlinear, t, u)
    T_tr = MM2015.T_triple(m.problem.thermo, eltype(u))
    return u[3] < T_tr || (u[3] == T_tr && non_latent_tendency(m, t, u) < 0)
end

function triple_value(m::Nonlinear, below, t, u)
    T_tr = MM2015.T_triple(m.problem.thermo, eltype(u))
    return below ? u[3] - T_tr : T_tr - u[3]
end

"""Whether the supersaturation of one phase (`liquid` or ice) rises with that phase inactive."""
function rising(m::FrozenLinear, liquid::Bool, other_active::Bool, t, u)
    a = liquid ? Active(false, other_active) : Active(other_active, false)
    return rhs(m, a, t, u)[1] > 0
end

function rising(m::Nonlinear, liquid::Bool, other_active::Bool, t, u)
    a = liquid ? Active(false, other_active) : Active(other_active, false)
    h = cbrt(eps(eltype(u)))
    F = rhs(m, a, t, u)
    index = liquid ? 1 : 2
    return supersaturations(m, t + h, u .+ h .* F)[index] - supersaturations(m, t - h, u .- h .* F)[index] > 0
end

"""
    activation(m, t, u) -> (Active, below)

Active phases at `(t, u)`: liquid when `x_l > 0` or `δ > 0`; below the triple point, ice when
`x_i > 0` or `δ_i > 0`; above it, ice when `x_i > 0` and `δ_i < 0`. A phase with zero mass exactly on
its saturation boundary is active when its supersaturation rises with the phase inactive.
"""
function activation(m, t, u)
    x_l, x_i = condensate(m, u)
    δ, δ_i = supersaturations(m, t, u)
    below = below_triple(m, t, u)
    liquid_plain = x_l > 0 || δ > 0
    ice_plain = below ? (x_i > 0 || δ_i > 0) : (x_i > 0 && δ_i < 0)
    liquid = liquid_plain || (iszero(x_l) && iszero(δ) && rising(m, true, ice_plain, t, u))
    ice = ice_plain || (below && iszero(x_i) && iszero(δ_i) && rising(m, false, liquid, t, u))
    return Active(liquid, ice), below
end

const EVENT_KINDS = (:liquid, :ice_mass, :ice_saturation, :triple)

"""
Indicator values of a segment with active phases `a`; the segment ends where one becomes positive:
`-x_l` or `δ` for liquid, `-x_i` and `δ_i` (or `-δ_i` for inactive ice with mass above the triple
point) for ice, and the signed distance past `T_triple`.
"""
function indicators(m, a::Active, below::Bool, t, u)
    x_l, x_i = condensate(m, u)
    δ, δ_i = supersaturations(m, t, u)
    none = -eltype(u)(Inf)
    liquid = a.liquid ? -x_l : δ
    ice_mass = a.ice ? -x_i : none
    ice_saturation = below ? (a.ice ? none : δ_i) : (a.ice ? δ_i : (x_i > 0 ? -δ_i : none))
    return (liquid, ice_mass, ice_saturation, triple_value(m, below, t, u))
end

fired(m, a, below, t, u) = any(>(0), indicators(m, a, below, t, u))

"""Radau IIA order-5 nodes and matrix in `FT`."""
function radau_tableau(::Type{FT}) where {FT}
    s6 = sqrt(FT(6))
    c = ((4 - s6) / 10, (4 + s6) / 10, one(FT))
    A = [
        (88 - 7s6)/360 (296 - 169s6)/1800 (-2 + 3s6)/225
        (296 + 169s6)/1800 (88 + 7s6)/360 (-2 - 3s6)/225
        (16 - s6)/36 (16 + s6)/36 one(FT)/9
    ]
    return c, A
end

"""In-place LU factorization with partial pivoting; returns the row permutation."""
function lu_factor!(M::Matrix)
    n = size(M, 1)
    perm = collect(1:n)
    for k in 1:n
        p = k - 1 + argmax(abs.(M[k:n, k]))
        if p != k
            M[[k, p], :] = M[[p, k], :]
            perm[[k, p]] = perm[[p, k]]
        end
        for i in (k + 1):n
            M[i, k] /= M[k, k]
            M[i, (k + 1):n] .-= M[i, k] .* M[k, (k + 1):n]
        end
    end
    return perm
end

function lu_solve(M::Matrix, perm, b::Vector)
    n = length(b)
    x = b[perm]
    for i in 2:n, j in 1:(i - 1)
        x[i] -= M[i, j] * x[j]
    end
    for i in n:-1:1
        for j in (i + 1):n
            x[i] -= M[i, j] * x[j]
        end
        x[i] /= M[i, i]
    end
    return x
end

"""Central-difference Jacobian of `rhs` with respect to the state."""
function jacobian(m, a, t, u, fd_scale)
    FT = eltype(u)
    n = length(u)
    J = zeros(FT, n, n)
    for j in 1:n
        h = cbrt(eps(FT)) * max(abs(u[j]), FT(fd_scale[j]))
        up = copy(u)
        um = copy(u)
        up[j] += h
        um[j] -= h
        J[:, j] = (rhs(m, a, t, up) - rhs(m, a, t, um)) / (2h)
    end
    return J
end

"""
One Radau IIA step of length `h` from `(t, u)` with the active phases `a`, solving the collocation
equations by simplified Newton iteration with the Jacobian `J`. The iteration stops when the update
is below `10⁻³` of the tolerance scale `atol + rtol |u|`, or below that scale and no longer halving
(the rounding floor of the right side). Returns `nothing` when neither happens in 100 iterations.
"""
function radau_step(m, a, t, u, h, J, tableau, atol, rtol)
    c, A = tableau
    FT = eltype(u)
    n = length(u)
    M = zeros(FT, 3n, 3n)
    for i in 1:3, j in 1:3, r in 1:n, s in 1:n
        M[(i - 1) * n + r, (j - 1) * n + s] = (i == j && r == s ? one(FT) : zero(FT)) - h * A[i, j] * J[r, s]
    end
    perm = lu_factor!(M)
    Z = zeros(FT, 3n)
    previous = FT(Inf)
    for _ in 1:100
        F = reduce(vcat, (rhs(m, a, t + c[i] * h, u .+ Z[((i - 1) * n + 1):(i * n)]) for i in 1:3))
        G = copy(Z)
        for i in 1:3, j in 1:3
            G[((i - 1) * n + 1):(i * n)] .-= h * A[i, j] .* F[((j - 1) * n + 1):(j * n)]
        end
        ΔZ = lu_solve(M, perm, -G)
        for i in 1:3, (active, liquid) in ((a.liquid, true), (a.ice, false))
            active || (ΔZ[(i - 1) * n + condensate_index(m, liquid)] = zero(FT))
        end
        Z .+= ΔZ
        update = maximum(1:(3n)) do k
            r = mod1(k, n)
            abs(ΔZ[k]) / (FT(atol[r]) + FT(rtol) * abs(u[r]))
        end
        (update ≤ FT(1e-3) || (update ≤ 1 && update ≥ previous / 2)) && return u .+ Z[(2n + 1):(3n)]
        previous = update
    end
    return nothing
end

"""Two Radau steps of length `h/2`: the step map used for accepted steps and event location."""
function double_step(m, a, t, u, h, J, tableau, atol, rtol)
    half = radau_step(m, a, t, u, h / 2, J, tableau, atol, rtol)
    half === nothing && return nothing, nothing
    return half, radau_step(m, a, t + h / 2, half, h / 2, J, tableau, atol, rtol)
end

"""
    Solution

Times `t` and states `u` at the ends of accepted steps and at events, the events as
`(kind, time)`, and the number of accepted steps.
"""
struct Solution{FT}
    t::Vector{FT}
    u::Vector{Vector{FT}}
    events::Vector{Tuple{Symbol, FT}}
    steps::Int
end

"""
    solve(model, Δt; rtol, atol_q, atol_T, h_max = Δt / 32, fixed_steps = 0, max_events = 1000) -> Solution

Integrate `model` over `[0, Δt]`. With `fixed_steps > 0`, take that many equal Radau steps and
require that no event occurs. With `fixed_steps = 0`, control the step by step doubling to `rtol`
and the absolute tolerances `atol_q` [kg kg⁻¹] and `atol_T` [K], and stop at every event.
"""
function solve(m, Δt::FT; rtol, atol_q, atol_T, h_max = Δt / 32, fixed_steps::Int = 0, max_events::Int = 1000) where {FT}
    atol, fd_scale = component_scales(m, atol_q, atol_T)
    tableau = radau_tableau(FT)
    t = zero(FT)
    u = initial_state(m)
    a, below = activation(m, t, u)
    ts, us, events = [t], [copy(u)], Tuple{Symbol, FT}[]
    steps = 0
    if fixed_steps > 0
        h = Δt / fixed_steps
        for _ in 1:fixed_steps
            J = jacobian(m, a, t, u, fd_scale)
            u_next = radau_step(m, a, t, u, h, J, tableau, atol, rtol)
            u_next === nothing && error("reference: Newton iteration failed at t = $t")
            t += h
            u = u_next
            steps += 1
            fired(m, a, below, t, u) && error("reference: event inside a fixed-step run at t = $t")
            push!(ts, t)
            push!(us, copy(u))
        end
        return Solution(ts, us, events, steps)
    end
    h = min(FT(h_max), Δt)
    while t < Δt
        h = min(h, Δt - t, FT(h_max))
        J = jacobian(m, a, t, u, fd_scale)
        full = radau_step(m, a, t, u, h, J, tableau, atol, rtol)
        half, doubled = double_step(m, a, t, u, h, J, tableau, atol, rtol)
        if full === nothing || doubled === nothing
            h /= 2
            continue
        end
        err = maximum(eachindex(u)) do r
            abs(doubled[r] - full[r]) / (31 * (FT(atol[r]) + FT(rtol) * max(abs(u[r]), abs(doubled[r]))))
        end
        if err > 1
            h *= max(FT(0.1), FT(0.9) * err^(-FT(1) / 6))
            continue
        end
        if fired(m, a, below, t + h, doubled)
            s_lo, s_hi = zero(FT), h
            while s_hi - s_lo > 4 * eps(FT) * (t + s_hi)
                s = (s_lo + s_hi) / 2
                (s == s_lo || s == s_hi) && break
                _, u_s = double_step(m, a, t, u, s, J, tableau, atol, rtol)
                u_s === nothing && error("reference: Newton iteration failed in event location at t = $(t + s)")
                fired(m, a, below, t + s, u_s) ? (s_hi = s) : (s_lo = s)
            end
            _, u_event = double_step(m, a, t, u, s_hi, J, tableau, atol, rtol)
            values = indicators(m, a, below, t + s_hi, u_event)
            t += s_hi
            u = u_event
            for (kind, value) in zip(EVENT_KINDS, values)
                value > 0 || continue
                push!(events, (kind, t))
                kind === :liquid && a.liquid && (u[condensate_index(m, true)] = zero(FT))
                kind === :ice_mass && (u[condensate_index(m, false)] = zero(FT))
            end
            length(events) > max_events && error("reference: more than $max_events events")
            a, below = activation(m, t, u)
        elseif fired(m, a, below, t + h / 2, half)
            h /= 2
            continue
        else
            t += h
            u = doubled
        end
        steps += 1
        push!(ts, t)
        push!(us, copy(u))
        h *= min(FT(4), FT(0.9) * max(err, eps(FT))^(-FT(1) / 6))
    end
    ts[end] = Δt
    return Solution(ts, us, events, steps)
end

"""Mean phase-change rates `(S̄_l, S̄_i)` over `[0, Δt]` in the basis of the model."""
function rates(m, sol::Solution, Δt)
    x_l0, x_i0 = condensate(m, first(sol.u))
    x_l1, x_i1 = condensate(m, last(sol.u))
    return ((x_l1 - x_l0) / Δt, (x_i1 - x_i0) / Δt)
end

end
