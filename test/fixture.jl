#=
Shared fixture for MM2015PiecewiseLinear / MM2015FixedT / MM2015 tests.
Fixture for Thermodynamics.jl
=#

using Test: Test
using ClimaParams
using Thermodynamics: Thermodynamics as TD
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

const FT = Float64

function mm2015_make_thermo_params(::Type{FT} = Float64) where {FT}
    return TD.Parameters.ThermodynamicsParameters(FT)
end

const MM2015_THERMO = mm2015_make_thermo_params()
const MM2015_T_FREEZE = FT(TD.Parameters.T_triple(MM2015_THERMO))

function mm2015_build_state(;
    freezing::Symbol = :BF,
    regime::Symbol = :wbf,
    q_liq::FT = FT(1e-6),
    q_ice::FT = FT(1e-6),
    τ_liq::FT = FT(10),
    τ_ice::FT = FT(10),
    Δt::FT = FT(10),
    w::FT = FT(0),
    dqvdt::FT = FT(0),
    dTdt::FT = FT(0),
    p::FT = FT(80000),
    sat_eps::FT = FT(1e-6),
    T::Union{FT, Nothing} = nothing,
    ρ::Union{FT, Nothing} = nothing,
    q_vap::Union{FT, Nothing} = nothing,
)
    T = if T !== nothing
        T
    elseif freezing === :BF
        FT(261.0)
    elseif freezing === :AF
        FT(280.0)
    elseif freezing === :near
        MM2015_T_FREEZE
    else
        error("unknown freezing=$freezing")
    end

    R_d = FT(TD.Parameters.R_d(MM2015_THERMO))
    ρ_seed = ρ !== nothing ? ρ : p / (R_d * T)
    q_vap_given = q_vap !== nothing

    if !q_vap_given
        _, _, q_vap = _mm2015_sat_and_qvapor(T, ρ_seed, regime, sat_eps)
    end
    q_tot = q_vap + q_liq + q_ice
    ρ = TD.air_density(MM2015_THERMO, T, p, q_tot, q_liq, q_ice)
    if !q_vap_given
        _, _, q_vap = _mm2015_sat_and_qvapor(T, ρ, regime, sat_eps)
        q_tot = q_vap + q_liq + q_ice
        ρ = TD.air_density(MM2015_THERMO, T, p, q_tot, q_liq, q_ice)
    end
    q_sl = TD.q_vap_saturation(MM2015_THERMO, T, ρ, TD.Liquid())
    q_si = TD.q_vap_saturation(MM2015_THERMO, T, ρ, TD.Ice())
    sc = MM2015.thermo_scalars(MM2015_THERMO, T, p, q_tot, q_liq, q_ice)
    return (;
        ρ,
        p,
        T,
        w,
        τ_liq,
        τ_ice,
        q_vap,
        q_tot,
        q_liq,
        q_ice,
        q_sl,
        q_si,
        Δt,
        dqvdt,
        dTdt,
        δ = q_vap - q_sl,
        δi = q_vap - q_si,
        thermo = MM2015_THERMO,
        sc,
    )
end

function _mm2015_sat_and_qvapor(T::FT, ρ::FT, regime::Symbol, sat_eps::FT)
    q_sl = TD.q_vap_saturation(MM2015_THERMO, T, ρ, TD.Liquid())
    q_si = TD.q_vap_saturation(MM2015_THERMO, T, ρ, TD.Ice())
    gap = q_sl - q_si
    q_vap = if regime === :super
        T < MM2015_T_FREEZE ? q_sl + abs(gap) * FT(0.5) + sat_eps : max(q_sl, q_si) + sat_eps
    elseif regime === :wbf
        (q_sl + q_si) / FT(2)
    elseif regime === :sub
        T < MM2015_T_FREEZE ? q_si - sat_eps : min(q_sl, q_si) - sat_eps
    elseif regime === :δ0
        q_sl
    elseif regime === :δi0
        q_si
    elseif regime === :δ_eps
        q_sl + eps(FT) * sign(gap == 0 ? one(FT) : gap)
    elseif regime === :δi_eps
        q_si + eps(FT)
    else
        error("unknown regime=$regime")
    end
    return q_sl, q_si, max(q_vap, FT(0))
end

function _mm2015_push!(states, kw)
    push!(states, mm2015_build_state(; kw...))
    return nothing
end

"""Designed state set by behavior class (not a blind Cartesian sweep)."""
function mm2015_designed_states(; include_extreme_τ::Bool = true)
    states = NamedTuple[]
    freezings = (:BF, :AF, :near)
    regimes_core = (:super, :wbf, :sub, :δ0, :δi0)
    presence_core = (
        (FT(0), FT(0)),
        (FT(1e-10), FT(0)),
        (FT(0), FT(1e-10)),
        (FT(1e-6), FT(0)),
        (FT(0), FT(1e-6)),
        (FT(1e-6), FT(1e-6)),
        (FT(1e-4), FT(1e-4)),
    )
    Δts_core = (FT(0.1), FT(10))
    τ_core = (
        (FT(1), FT(1)),
        (FT(0.01), FT(0.01)),
        (FT(100), FT(100)),
        (FT(0.01), FT(1000)),
        (FT(1000), FT(0.01)),
    )
    τ_extreme = ((FT(1e-12), FT(1e9)), (FT(1e9), FT(1e-12)), (FT(1), FT(5e5)), (FT(80), FT(5e5)))

    for fr in freezings, rg in regimes_core, (ql, qi) in presence_core, Δt in Δts_core, (τl, τi) in τ_core
        _mm2015_push!(states, (; freezing = fr, regime = rg, q_liq = ql, q_ice = qi, τ_liq = τl, τ_ice = τi, Δt = Δt))
    end
    for fr in (:BF, :AF), rg in (:δ_eps, :δi_eps)
        _mm2015_push!(states, (; freezing = fr, regime = rg, q_liq = FT(1e-6), q_ice = FT(1e-6), τ_liq = FT(1), τ_ice = FT(1)))
    end
    if include_extreme_τ
        for fr in (:BF, :AF), rg in (:sub, :wbf), (ql, qi) in ((FT(1e-14), FT(1e-14)), (FT(1e-10), FT(4e-11))), (τl, τi) in τ_extreme
            _mm2015_push!(states, (; freezing = fr, regime = rg, q_liq = ql, q_ice = qi, τ_liq = τl, τ_ice = τi, Δt = FT(10), sat_eps = FT(1e-4)))
        end
    end
    for fr in (:BF, :AF), rg in (:sub, :wbf)
        _mm2015_push!(states, (; freezing = fr, regime = rg, q_liq = FT(1e-6), q_ice = FT(1e-6), τ_liq = FT(10), τ_ice = FT(100), w = FT(0.1)))
        _mm2015_push!(states, (; freezing = fr, regime = rg, q_liq = FT(1e-6), q_ice = FT(1e-6), τ_liq = FT(10), τ_ice = FT(100), dqvdt = FT(1e-6), dTdt = FT(1e-4)))
    end
    return states
end

function mm2015_inputs(st)
    sc = st.sc
    return MM2015.MM2015Inputs(
        sc.g, sc.L_i, sc.L_l, sc.c_p, sc.T_freeze,
        sc.dqsl_dT, sc.dqsi_dT, sc.e_sl, sc.e_si,
        st.ρ, st.p, st.T, st.w, st.τ_liq, st.τ_ice,
        st.q_tot, st.q_liq, st.q_ice, st.q_sl, st.q_si, st.Δt,
    )
end

function mm2015_call_sources(scheme, st; opts = nothing)
    inputs = mm2015_inputs(st)
    if scheme isa MM2015.MM2015PiecewiseLinear
        o = opts === nothing ? MM2015.MM2015PiecewiseLinearOpts{FT}() : opts
        return MM2015.morrison_milbrandt_2015(scheme, inputs; opts = o, dqvdt = st.dqvdt, dTdt = st.dTdt)
    elseif scheme isa MM2015.MM2015FixedT
        o = opts === nothing ? MM2015.MM2015FixedTOpts{FT}(; dqvdt = st.dqvdt, dTdt = st.dTdt) : opts
        return MM2015.morrison_milbrandt_2015(scheme, inputs; opts = o)
    elseif scheme isa MM2015.MM2015
        o = opts === nothing ? MM2015.MM2015Opts{FT}(; dqvdt = st.dqvdt, dTdt = st.dTdt) : opts
        return MM2015.morrison_milbrandt_2015(scheme, inputs; opts = o, thermo = st.thermo)
    else
        error("unknown scheme $scheme")
    end
end

function mm2015_call_epa_solver(st; opts = MM2015.MM2015FixedTOpts{FT}(; dqvdt = st.dqvdt, dTdt = st.dTdt))
    return mm2015_call_sources(MM2015.MM2015FixedT(), st; opts = opts)
end

function mm2015_assert_invariants(S_ql, S_qi, st; atol = eps(FT))
    Test.@test isfinite(S_ql) && isfinite(S_qi)
    ql, qi = st.q_liq, st.q_ice
    Δt = st.Δt
    iszero(ql) && Test.@test S_ql ≥ -atol
    iszero(qi) && Test.@test S_qi ≥ -atol
    st.T > MM2015_T_FREEZE && Test.@test S_qi ≤ atol
    Test.@test S_ql * Δt ≥ -(ql + atol)
    Test.@test S_qi * Δt ≥ -(qi + atol)
    Test.@test (S_ql + S_qi) * Δt ≥ -(ql + qi + atol)
    q_vap_avail = st.q_vap + max(st.dqvdt, zero(FT)) * Δt + atol
    Test.@test (S_ql + S_qi) * Δt ≤ q_vap_avail + FT(1e-10)
end
