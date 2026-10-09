# Designed parcel problems shared by tests, benchmarks, and figures.

module ParcelCorpus

using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

"""
Corpus cases in specific humidity: temperature `T` [K], pressure `p` [Pa], liquid `q_l` and ice
`q_i` [kg kg⁻¹], vapor set by `humidity`, relaxation times `τ_l` and `τ_i` [s], vertical velocity
`w` [m s⁻¹] (`:band` for `band_velocity`), temperature forcing `dTdt` [K s⁻¹], vapor
forcing `dq_vap_dt` [s⁻¹], and step `Δt` [s]. `humidity` is `(:RH_l, x)`, `(:RH_i, x)`,
`(:δ, x)`, `(:δ_i, x)`, or `(:wbf, x)` for `δ = -x Δ`.
"""
const CORPUS = (
    (; name = :warm_updraft, T = 285.0, p = 9.0e4, humidity = (:RH_l, 1.001), q_l = 2e-4, q_i = 0.0, τ_l = 5.0, τ_i = 1e3, w = 1.0, dTdt = 0.0, dq_vap_dt = 0.0, Δt = 60.0),
    (; name = :warm_evaporation, T = 290.0, p = 9.5e4, humidity = (:RH_l, 0.9), q_l = 1e-5, q_i = 0.0, τ_l = 10.0, τ_i = 1e3, w = 0.0, dTdt = 0.0, dq_vap_dt = 0.0, Δt = 120.0),
    (; name = :wbf, T = 261.0, p = 8.0e4, humidity = (:RH_l, 0.995), q_l = 2e-4, q_i = 1e-4, τ_l = 8.0, τ_i = 60.0, w = 0.5, dTdt = 0.0, dq_vap_dt = 0.0, Δt = 300.0),
    (; name = :ice_subliming_liquid_growing, T = 275.0, p = 8.5e4, humidity = (:RH_l, 1.0005), q_l = 1e-4, q_i = 5e-5, τ_l = 10.0, τ_i = 20.0, w = 0.0, dTdt = 0.0, dq_vap_dt = 0.0, Δt = 120.0),
    (; name = :ice_only_supersaturated, T = 225.0, p = 2.5e4, humidity = (:RH_i, 1.2), q_l = 0.0, q_i = 1e-5, τ_l = 1e3, τ_i = 600.0, w = 0.2, dTdt = 0.0, dq_vap_dt = 0.0, Δt = 600.0),
    (; name = :activation_moistening, T = 280.0, p = 9.0e4, humidity = (:RH_l, 0.999), q_l = 0.0, q_i = 0.0, τ_l = 5.0, τ_i = 1e3, w = 0.0, dTdt = 0.0, dq_vap_dt = 2e-6, Δt = 30.0),
    (; name = :activation_ascent, T = 280.0, p = 9.0e4, humidity = (:RH_l, 0.999), q_l = 0.0, q_i = 0.0, τ_l = 5.0, τ_i = 1e3, w = 1.0, dTdt = 0.0, dq_vap_dt = 0.0, Δt = 30.0),
    (; name = :activation_cooling, T = 280.0, p = 9.0e4, humidity = (:RH_l, 0.999), q_l = 0.0, q_i = 0.0, τ_l = 5.0, τ_i = 1e3, w = 0.0, dTdt = -1e-3, dq_vap_dt = 0.0, Δt = 30.0),
    (; name = :freezing_level, T = 273.4, p = 7.0e4, humidity = (:RH_l, 1.0), q_l = 1e-4, q_i = 1e-5, τ_l = 10.0, τ_i = 100.0, w = 0.0, dTdt = -2e-3, dq_vap_dt = 0.0, Δt = 300.0),
    (; name = :stiff, T = 280.0, p = 9.0e4, humidity = (:RH_l, 1.002), q_l = 5e-4, q_i = 0.0, τ_l = 0.1, τ_i = 1e3, w = 2.0, dTdt = 0.0, dq_vap_dt = 0.0, Δt = 60.0),
    (; name = :sluggish, T = 280.0, p = 9.0e4, humidity = (:RH_l, 1.01), q_l = 1e-6, q_i = 0.0, τ_l = 1e4, τ_i = 1e4, w = 0.0, dTdt = 0.0, dq_vap_dt = 0.0, Δt = 10.0),
    (; name = :dual_depletion, T = 265.0, p = 8.0e4, humidity = (:RH_l, 0.8), q_l = 1e-6, q_i = 1e-6, τ_l = 5.0, τ_i = 20.0, w = 0.0, dTdt = 0.0, dq_vap_dt = 0.0, Δt = 120.0),
    (; name = :band, T = 261.0, p = 8.0e4, humidity = (:δ, -1e-7), q_l = 0.0, q_i = 1e-4, τ_l = 10.0, τ_i = 50.0, w = :band, dTdt = 0.0, dq_vap_dt = 0.0, Δt = 60.0),
    (; name = :massless_ice_wbf, T = 261.81, p = 8.0e4, humidity = (:wbf, 0.5), q_l = 2e-6, q_i = 0.0, τ_l = 10.0, τ_i = 80.0, w = 0.0, dTdt = 0.0, dq_vap_dt = 0.0, Δt = 10.0),
    (; name = :ice_subliming_above_freezing, T = 276.77, p = 9.0e4, humidity = (:RH_l, 0.98), q_l = 0.0, q_i = 7e-4, τ_l = 6000.0, τ_i = 4.0, w = 0.0, dTdt = 0.0, dq_vap_dt = 0.0, Δt = 10.0),
    (; name = :float_scale_near_freeze, T = 269.4, p = 88268.0, humidity = (:δ, -1.6e-5), q_l = 1e-10, q_i = 4e-11, τ_l = 80.0, τ_i = 5e5, w = 0.0, dTdt = 2.8e-4, dq_vap_dt = 3.6e-6, Δt = 10.0),
    (; name = :float_scale_dual_depletion, T = 270.2, p = 88268.0, humidity = (:δ_i, -2.384615384615381e-5), q_l = 1e-10, q_i = 1e-10, τ_l = 80.0, τ_i = 5e5, w = 0.0, dTdt = 2.8e-4, dq_vap_dt = 3.6e-6, Δt = 10.0),
)

"""The corpus case named `name`."""
corpus_case(name::Symbol) = only(filter(c -> c.name == name, CORPUS))

"""
    corpus_vapor(thermo, T, p, q_l, q_i, humidity) -> q_v

Vapor specific humidity that meets `humidity` with pressure-based saturation `q_s = (1 − q_t) r_s`.
"""
function corpus_vapor(thermo, T::FT, p::FT, q_l::FT, q_i::FT, humidity) where {FT}
    kind, value = humidity
    r_l = MM2015.saturation(thermo, T, p, MM2015.Liquid()).r
    r_i = MM2015.saturation(thermo, T, p, MM2015.Ice()).r
    q_c = q_l + q_i
    x = FT(value)
    kind === :RH_l && return x * r_l * (1 - q_c) / (1 + x * r_l)
    kind === :RH_i && return x * r_i * (1 - q_c) / (1 + x * r_i)
    kind === :δ && return (x + r_l * (1 - q_c)) / (1 + r_l)
    kind === :δ_i && return (x + r_i * (1 - q_c)) / (1 + r_i)
    kind === :wbf && return ((1 - x) * r_l + x * r_i) * (1 - q_c) / (1 + (1 - x) * r_l + x * r_i)
    throw(ArgumentError("unknown humidity kind $kind"))
end

"""
    band_velocity(thermo, T, p, q_t, q_l, q_i, τ_i) -> w

Vertical velocity at the centre of the band where, with liquid absent and forcing by ascent only,
the ice-frame forcing `A_i` lies in `(Δ/τ_i, αΔ/(τ_iΓ_i))` and `0 < A_l < αΔ/(τ_iΓ_i)`.
"""
function band_velocity(thermo, T::FT, p::FT, q_t::FT, q_l::FT, q_i::FT, τ_i::FT) where {FT}
    unit = MM2015.MM2015Problem(
        MM2015.SpecificHumidity(), thermo, MM2015.MM2015State(T, p, q_t, q_l, q_i),
        MM2015.MM2015Timescales(FT(1), τ_i),
        MM2015.MM2015Forcing(-MM2015.air_density(thermo, T, p, q_t, q_l, q_i) * MM2015.grav(thermo, FT), zero(FT), zero(FT)),
    )
    k = MM2015.coefficients(unit)
    ice = MM2015.saturation(thermo, T, p, MM2015.Ice())
    q_d = 1 - q_t
    ρ = MM2015.air_density(thermo, T, p, q_t, q_l, q_i)
    c_pm = MM2015.cp_m(thermo, q_t, q_l, q_i)
    dpdt = unit.forcing.dpdt
    A_i = -q_d * ice.dr_dT * dpdt / (ρ * c_pm) - q_d * ice.dr_dp * dpdt
    upper = k.α * k.Δ / (τ_i * k.Γ_i)
    lo = k.Δ / (τ_i * A_i)
    hi = min(upper / A_i, upper / k.A_l)
    lo < hi || throw(ArgumentError("empty band at T = $T, p = $p"))
    return (lo + hi) / 2
end

"""
    corpus_problem(case, thermo; basis = SpecificHumidity(), FT = Float64) -> (problem, Δt)

The `MM2015Problem` of a corpus `case` in `basis` and floating-point type `FT`.
"""
function corpus_problem(case, thermo; basis = MM2015.SpecificHumidity(), FT = Float64)
    T, p, q_l, q_i = FT(case.T), FT(case.p), FT(case.q_l), FT(case.q_i)
    q_v = corpus_vapor(thermo, T, p, q_l, q_i, case.humidity)
    q_t = q_v + q_l + q_i
    w = case.w === :band ? band_velocity(thermo, T, p, q_t, q_l, q_i, FT(case.τ_i)) : FT(case.w)
    dpdt = -MM2015.air_density(thermo, T, p, q_t, q_l, q_i) * MM2015.grav(thermo, FT) * w
    F_v = FT(case.dq_vap_dt)
    state, dx_vap_dt = if basis isa MM2015.SpecificHumidity
        MM2015.MM2015State(T, p, q_t, q_l, q_i), F_v
    else
        q_d = 1 - q_t
        MM2015.MM2015State(T, p, q_t / q_d, q_l / q_d, q_i / q_d), F_v / q_d + q_v * F_v / q_d^2
    end
    problem = MM2015.MM2015Problem(
        basis, thermo, state,
        MM2015.MM2015Timescales(FT(case.τ_l), FT(case.τ_i)),
        MM2015.MM2015Forcing(dpdt, FT(case.dTdt), dx_vap_dt),
    )
    return problem, FT(case.Δt)
end

end
