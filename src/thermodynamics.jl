#=

Thermodynamics backends for residual T-updating ([`MM2015`](@ref)).

The pipeline uses scalar fields `(T, p, q_tot, q_liq, q_ice)`. There is no stored
thermodynamic state. Dispatch is on the `thermo` handle:

- [`DefaultThermodynamicsBackend`](@ref) — ClimaParams `parameters.toml` numbers,
  Clausius–Clapeyron with constant `L_v0`/`L_s0`, Tripoli–Cotton `θ_li`.
- Thermodynamics.jl via MorrisonMilbrandt2015ThermodynamicsExt (`params` plus variables).

=#

"""
    AbstractThermodynamicsBackend

Supertype for the built-in, dependency-free thermodynamics backend. Extension backends dispatch
on their own parameter set and need not subtype this.
"""
abstract type AbstractThermodynamicsBackend end
Base.broadcastable(backend::AbstractThermodynamicsBackend) = tuple(backend)

"""
    DefaultThermodynamicsBackend()

Built-in backend: ideal-gas density, Clausius–Clapeyron saturation (constant `L_v0`/`L_s0`),
Tripoli–Cotton `θ_li`. Constants match ClimaParams.jl `src/parameters.toml` (no ClimaParams
dependency).
"""
struct DefaultThermodynamicsBackend <: AbstractThermodynamicsBackend end

"""Moist-air isobaric heat capacity [J/kg/K], same `c_p` as `θ_li`."""
function cp_m end
"""Moist air density [kg/m³] at prescribed condensate (noneq)."""
function air_density end
"""Virtual temperature [K] at prescribed condensate."""
function virtual_temperature end
"""Liquid-ice potential temperature [K] at prescribed condensate."""
function liquid_ice_pottemp end
"""Saturation vapor pressure [Pa] for the given phase."""
function saturation_vapor_pressure end
"""Saturation vapor pressure over liquid water [Pa]."""
function saturation_vapor_pressure_liq end
"""Saturation vapor pressure over ice [Pa]."""
function saturation_vapor_pressure_ice end
"""Dry-air gas constant [J/kg/K]."""
function R_d end
"""Water-vapor gas constant [J/kg/K]."""
function R_v end
"""Gravitational acceleration [m/s^2]."""
function grav end
"""Dry-air isobaric heat capacity [J/kg/K]."""
function cp_d end
"""Water-vapor isobaric heat capacity [J/kg/K]."""
function cp_v end
"""Liquid-water isobaric heat capacity [J/kg/K]."""
function cp_l end
"""Ice isobaric heat capacity [J/kg/K]."""
function cp_i end
"""Molar mass of dry air [kg/mol]."""
function molar_mass_dry_air end
"""Molar mass of water [kg/mol]."""
function molar_mass_water end
"""Ratio of molar masses Mv/Md."""
function molmass_ratio end
"""Reference pressure for potential temperature [Pa]."""
function p_ref end
"""Freezing temperature [K] (melting point). Freeze *events* use [`T_triple`](@ref)."""
function T_freeze end
"""Triple-point temperature [K]. MM2015 freeze events use this."""
function T_triple end
"""Latent heat of vaporization at the reference temperature [J/kg]."""
function L_v0 end
"""Latent heat of sublimation at the reference temperature [J/kg]."""
function L_s0 end
"""Saturation vapor pressure at the triple point [Pa]; the CC anchor for both phases."""
function e_ref end
"""Reference temperature the latent heats are anchored at [K] (`T_0` = `T_triple`)."""
function T_0 end
"""Latent heat at `T` [J/kg] from its reference value `LH_0` and the heat-capacity difference `Δcp`."""
function latent_heat_generic end
"""Latent heat of vaporization at `T` [J/kg]."""
function latent_heat_vapor end
"""Latent heat of sublimation at `T` [J/kg]."""
function latent_heat_sublim end
"""∂q_vap_sat_liq/∂T [kg/kg/K]."""
function ∂q_vap_sat_liq_∂T end
"""∂q_vap_sat_ice/∂T [kg/kg/K]."""
function ∂q_vap_sat_ice_∂T end
"""NonEquilibrium temperature from `(p, θ_li)` at fixed `q_liq`, `q_ice`."""
function air_temperature_noneq_pθq end
"""`(∂θ_li/∂T)_{p,q}` used to evolve `θ_li` under host `dTdt`."""
function dθ_li_dT end
"""Latent heats, `c_p`, `∂q*/∂T`, saturation, and density at `(T, p, q)`."""
function thermo_scalars end
"""Moist-air gas constant `R_m = R_d (1 − q_tot) + R_v q_vap` [J/kg/K]."""
function gas_constant_air end
"""`dR_m/dt` along a trajectory with rates `q_tot'`, `q_liq'`, `q_ice'`."""
function d_gas_constant_air_dt end
"""`dT/dt` from the noneq `pθq` invert along a trajectory."""
function dT_noneq_dt end
"""`d q_sl / dt` along a trajectory (`q_sl` over liquid)."""
function dq_vap_sat_liq_dt end
"""`d q_si / dt` along a trajectory (`q_si` over ice)."""
function dq_vap_sat_ice_dt end

"""Supertype for the water phase a saturation quantity is taken over."""
abstract type AbstractPhase end
"""Water vapor."""
struct Vapor <: AbstractPhase end
"""Liquid water; saturation is taken over a plane liquid surface."""
struct Liquid <: AbstractPhase end
"""Ice; saturation is taken over a plane ice surface."""
struct Ice <: AbstractPhase end

# ClimaParams.jl src/parameters.toml (no ClimaParams dependency).
@inline e_ref(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(611.657)
@inline R_d(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(287.0)
@inline R_v(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(461.5)
@inline grav(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(9.81)
@inline cp_d(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(1004.5)
@inline cp_v(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(1859)
@inline cp_l(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(4181)
@inline cp_i(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(2070)
@inline molar_mass_dry_air(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(0.02897)
@inline molar_mass_water(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(0.01801528)
@inline molmass_ratio(backend::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} =
    molar_mass_water(backend, FT) / molar_mass_dry_air(backend, FT)
@inline p_ref(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(1.0e5)
@inline T_0(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(273.16)
@inline T_triple(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(273.16)
@inline T_freeze(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(273.15)
@inline L_v0(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(2.5008e6)
@inline L_s0(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(2.8344e6)

@inline function cp_m(
    backend::DefaultThermodynamicsBackend,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
) where {FT}
    return cp_d(backend, FT) +
           (cp_v(backend, FT) - cp_d(backend, FT)) * q_tot +
           (cp_l(backend, FT) - cp_v(backend, FT)) * q_liq +
           (cp_i(backend, FT) - cp_v(backend, FT)) * q_ice
end

@inline latent_heat_generic(backend::DefaultThermodynamicsBackend, T::FT, LH_0::FT, Δcp::FT) where {FT} =
    LH_0 + Δcp * (T - T_0(backend, FT))
@inline latent_heat_vapor(backend::DefaultThermodynamicsBackend, T::FT) where {FT} =
    latent_heat_generic(backend, T, L_v0(backend, FT), cp_v(backend, FT) - cp_l(backend, FT))
@inline latent_heat_sublim(backend::DefaultThermodynamicsBackend, T::FT) where {FT} =
    latent_heat_generic(backend, T, L_s0(backend, FT), cp_v(backend, FT) - cp_i(backend, FT))

"""
    saturation_vapor_pressure_liq(backend, T)

Saturation vapor pressure over liquid [Pa], Clausius–Clapeyron with constant `L_v0`,
anchored at `(T_triple, press_triple)`.
"""
@inline function saturation_vapor_pressure_liq(
    backend::DefaultThermodynamicsBackend,
    T::FT;
    T_fr::FT = T_triple(backend, FT),
    e_anchor::FT = e_ref(backend, FT),
) where {FT}
    return e_anchor * exp(L_v0(backend, FT) / R_v(backend, FT) * (one(FT) / T_fr - one(FT) / T))
end

"""
    saturation_vapor_pressure_ice(backend, T)

Saturation vapor pressure over ice [Pa], Clausius–Clapeyron with constant `L_s0`,
anchored at `(T_triple, press_triple)`.
"""
@inline function saturation_vapor_pressure_ice(
    backend::DefaultThermodynamicsBackend,
    T::FT;
    T_fr::FT = T_triple(backend, FT),
    e_anchor::FT = e_ref(backend, FT),
) where {FT}
    return e_anchor * exp(L_s0(backend, FT) / R_v(backend, FT) * (one(FT) / T_fr - one(FT) / T))
end

@inline function saturation_vapor_pressure(
    backend::DefaultThermodynamicsBackend,
    T::FT,
    phase::AbstractPhase = Liquid();
    T_fr::FT = T_triple(backend, FT),
    e_anchor::FT = e_ref(backend, FT),
) where {FT}
    if phase === Liquid()
        return saturation_vapor_pressure_liq(backend, T; T_fr, e_anchor)
    elseif phase === Ice()
        return saturation_vapor_pressure_ice(backend, T; T_fr, e_anchor)
    else
        error("saturation_vapor_pressure: phase must be Liquid() or Ice()")
    end
end

@inline vapor_specific_humidity(q_tot::FT, q_liq::FT, q_ice::FT) where {FT} = q_tot - q_liq - q_ice
@inline vapor_specific_humidity(::DefaultThermodynamicsBackend, q_tot::FT, q_liq::FT, q_ice::FT) where {FT} =
    vapor_specific_humidity(q_tot, q_liq, q_ice)

@inline function virtual_temperature(
    backend::DefaultThermodynamicsBackend,
    T::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT;
    R_v::FT = R_v(backend, FT),
    R_d::FT = R_d(backend, FT),
) where {FT}
    q_vap = q_tot - q_liq - q_ice
    return T * (one(FT) + (R_v / R_d - one(FT)) * q_vap - q_liq - q_ice)
end

@inline function air_density(
    backend::DefaultThermodynamicsBackend,
    T::FT,
    p::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT;
    R_v::FT = R_v(backend, FT),
    R_d::FT = R_d(backend, FT),
) where {FT}
    return p / (R_d * virtual_temperature(backend, T, q_tot, q_liq, q_ice; R_v, R_d))
end

"""
    liquid_ice_pottemp(backend, T, p, q_tot, q_liq, q_ice)

Tripoli–Cotton liquid-ice potential temperature at **prescribed** condensate:

`θ_li = T Π exp(−A/T)`, `Π = (p_ref/p)^{R_d/c_p}`, `A = (L_v0 q_liq + L_s0 q_ice)/c_p`,
with moist `c_p = cp_m`.
"""
@inline function liquid_ice_pottemp(
    backend::DefaultThermodynamicsBackend,
    T::FT,
    p::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
) where {FT}
    cpm = cp_m(backend, q_tot, q_liq, q_ice)
    κ = R_d(backend, FT) / cpm
    A = (L_v0(backend, FT) * q_liq + L_s0(backend, FT) * q_ice) / cpm
    return T * (p_ref(backend, FT) / p)^κ * exp(-A / T)
end

"""`(∂θ_li/∂T)_{p,q} = (θ/T)(1 + A/T)` for Tripoli–Cotton `θ_li`."""
@inline function dθ_li_dT(
    backend::DefaultThermodynamicsBackend,
    T::FT,
    p::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
) where {FT}
    θ = liquid_ice_pottemp(backend, T, p, q_tot, q_liq, q_ice)
    cpm = cp_m(backend, q_tot, q_liq, q_ice)
    A = (L_v0(backend, FT) * q_liq + L_s0(backend, FT) * q_ice) / cpm
    return (θ / T) * (one(FT) + A / T)
end

"""`θ_li(t) = θ_li(0) + (∂θ_li/∂T)_{p,q} dTdt t` on a frozen-coeff segment."""
@inline θ_li_evolved(θ_li0, dθdT0, dTdt, t) = θ_li0 + dθdT0 * dTdt * t

"""
    air_temperature_noneq_pθq(backend, p, θ_li, q_tot, q_liq, q_ice)

Invert Tripoli–Cotton `θ_li` at **fixed** condensate (NonEquilibrium). Lambert W₀;
not saturation adjustment. `A = 0` ⇒ `T = θ_li / Π`.
"""
function air_temperature_noneq_pθq(
    backend::DefaultThermodynamicsBackend,
    p::FT,
    θ_li::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
) where {FT}
    cpm = cp_m(backend, q_tot, q_liq, q_ice)
    κ = R_d(backend, FT) / cpm
    Π = (p_ref(backend, FT) / p)^κ
    A = (L_v0(backend, FT) * q_liq + L_s0(backend, FT) * q_ice) / cpm
    iszero(A) && return θ_li / Π
    return A / fast_lambertw0(A * Π / θ_li)
end

@inline function ∂q_from_e_p(e_sat::FT, p::FT, ε::FT, L::FT, T::FT, R_v::FT) where {FT}
    denom = p - (one(FT) - ε) * e_sat
    denom > zero(FT) || return zero(FT)
    dq_de = ε * p / (denom^2)
    de_dT = e_sat * L / (R_v * T^2)
    return dq_de * de_dT
end

function ∂q_vap_sat_liq_∂T(backend::DefaultThermodynamicsBackend, T::FT, p::FT) where {FT}
    e_sat = saturation_vapor_pressure_liq(backend, T)
    return ∂q_from_e_p(e_sat, p, molmass_ratio(backend, FT), L_v0(backend, FT), T, R_v(backend, FT))
end

function ∂q_vap_sat_ice_∂T(backend::DefaultThermodynamicsBackend, T::FT, p::FT) where {FT}
    e_sat = saturation_vapor_pressure_ice(backend, T)
    return ∂q_from_e_p(e_sat, p, molmass_ratio(backend, FT), L_s0(backend, FT), T, R_v(backend, FT))
end

@inline function thermo_scalars(
    backend::DefaultThermodynamicsBackend,
    T::FT,
    p::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
) where {FT}
    ρ = air_density(backend, T, p, q_tot, q_liq, q_ice)
    e_sl = saturation_vapor_pressure_liq(backend, T)
    e_si = saturation_vapor_pressure_ice(backend, T)
    ε = molmass_ratio(backend, FT)
    q_d = one(FT) - q_tot
    r_sl = ε * e_sl / (p - e_sl)
    r_si = ε * e_si / (p - e_si)
    q_sl = q_d * r_sl
    q_si = q_d * r_si
    drsl_dT = ε * p / (p - e_sl)^2 * e_sl * L_v0(backend, FT) / (R_v(backend, FT) * T^2)
    drsi_dT = ε * p / (p - e_si)^2 * e_si * L_s0(backend, FT) / (R_v(backend, FT) * T^2)
    return (
        g = grav(backend, FT),
        L_i = latent_heat_sublim(backend, T),
        L_l = latent_heat_vapor(backend, T),
        c_p = cp_m(backend, q_tot, q_liq, q_ice),
        T_freeze = T_triple(backend, FT),
        dqsl_dT = q_d * drsl_dT,
        dqsi_dT = q_d * drsi_dT,
        e_sl,
        e_si,
        ρ,
        q_sl,
        q_si,
        q_vap = vapor_specific_humidity(q_tot, q_liq, q_ice),
    )
end

@inline function gas_constant_air(
    backend::DefaultThermodynamicsBackend,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
) where {FT}
    q_vap = q_tot - q_liq - q_ice
    return R_d(backend, FT) * (one(FT) - q_tot) + R_v(backend, FT) * q_vap
end

@inline function d_gas_constant_air_dt(
    backend::DefaultThermodynamicsBackend,
    dqt::FT,
    dql::FT,
    dqi::FT,
) where {FT}
    return -R_d(backend, FT) * dqt + R_v(backend, FT) * (dqt - dql - dqi)
end

@inline function ∂q_vap_sat_liq_∂p(backend::DefaultThermodynamicsBackend, T::FT, p::FT) where {FT}
    e_sat = saturation_vapor_pressure_liq(backend, T)
    ε = molmass_ratio(backend, FT)
    denom = p - (one(FT) - ε) * e_sat
    denom > zero(FT) || return zero(FT)
    return -ε * e_sat / (denom^2)
end

@inline function ∂q_vap_sat_ice_∂p(backend::DefaultThermodynamicsBackend, T::FT, p::FT) where {FT}
    e_sat = saturation_vapor_pressure_ice(backend, T)
    ε = molmass_ratio(backend, FT)
    denom = p - (one(FT) - ε) * e_sat
    denom > zero(FT) || return zero(FT)
    return -ε * e_sat / (denom^2)
end

"""Derivative of paper-based `q_sl=q_d εe_sl/(p-e_sl)` along a parcel trajectory."""
function dq_vap_sat_liq_dt(
    backend::DefaultThermodynamicsBackend,
    T::FT,
    p::FT,
    q_tot::FT,
    q_sl::FT,
    dTdt::FT,
    dpdt::FT,
    dqtotdt::FT,
) where {FT}
    q_d = one(FT) - q_tot
    r_sl = q_sl / q_d
    e_sl = saturation_vapor_pressure_liq(backend, T)
    drdT = molmass_ratio(backend, FT) * p / (p - e_sl)^2 *
            e_sl * L_v0(backend, FT) / (R_v(backend, FT) * T^2)
    drdp = -r_sl / (p - e_sl)
    return -dqtotdt * r_sl + q_d * (drdT * dTdt + drdp * dpdt)
end

"""Derivative of paper-based `q_si=q_d εe_si/(p-e_si)` along a parcel trajectory."""
function dq_vap_sat_ice_dt(
    backend::DefaultThermodynamicsBackend,
    T::FT,
    p::FT,
    q_tot::FT,
    q_si::FT,
    dTdt::FT,
    dpdt::FT,
    dqtotdt::FT,
) where {FT}
    q_d = one(FT) - q_tot
    r_si = q_si / q_d
    e_si = saturation_vapor_pressure_ice(backend, T)
    drdT = molmass_ratio(backend, FT) * p / (p - e_si)^2 *
            e_si * L_s0(backend, FT) / (R_v(backend, FT) * T^2)
    drdp = -r_si / (p - e_si)
    return -dqtotdt * r_si + q_d * (drdT * dTdt + drdp * dpdt)
end

"""
    dT_noneq_dt(backend, T, p, θ_li, q_tot, q_liq, q_ice, dpdt, dθdt, dqt, dql, dqi)

Chain rule for the Tripoli–Cotton invert `T = A / W₀(A Π / θ)` (`A = 0` ⇒ `T = θ/Π`).

`Π = (p_ref/p)^{R_d/c_p}`, `A = (L_v0 q_liq + L_s0 q_ice)/c_p`. For `A ≠ 0`,
`W = A/T` and

`T' / T = (A'/A) W/(1+W) − (Π'/Π)/(1+W) + (θ'/θ)/(1+W)`.
"""
function dT_noneq_dt(
    backend::DefaultThermodynamicsBackend,
    T::FT,
    p::FT,
    θ_li::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
    dpdt::FT,
    dθdt::FT,
    dqt::FT,
    dql::FT,
    dqi::FT,
) where {FT}
    cpm = cp_m(backend, q_tot, q_liq, q_ice)
    dcpm =
        (cp_v(backend, FT) - cp_d(backend, FT)) * dqt +
        (cp_l(backend, FT) - cp_v(backend, FT)) * dql +
        (cp_i(backend, FT) - cp_v(backend, FT)) * dqi
    κ = R_d(backend, FT) / cpm
    dκ = -κ / cpm * dcpm
    ln_pr = log(p_ref(backend, FT) / p)
    dlnΠ = dκ * ln_pr - (κ / p) * dpdt
    A = (L_v0(backend, FT) * q_liq + L_s0(backend, FT) * q_ice) / cpm
    dA = (L_v0(backend, FT) * dql + L_s0(backend, FT) * dqi) / cpm - A / cpm * dcpm
    if iszero(A)
        return T * (dθdt / θ_li - dlnΠ) + dA
    end
    W = A / T
    return T * (dA / A * W + (-dlnΠ + dθdt / θ_li)) / (one(FT) + W)
end

@inline function dρ_dt(ρ::FT, p::FT, T::FT, R_m::FT, dpdt::FT, dTdt::FT, dRmdt::FT) where {FT}
    return ρ * (dpdt / p - dTdt / T - dRmdt / R_m)
end
