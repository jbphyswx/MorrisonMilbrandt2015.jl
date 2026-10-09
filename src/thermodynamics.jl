"""
    DefaultThermodynamicsBackend()

Thermodynamics backend without dependencies, with fixed gas constants, heat capacities, gravity, and
triple point. Its saturation vapor pressures integrate the Clausius–Clapeyron relation with latent heats that vary
linearly with temperature (Kirchhoff), anchored at the triple point.
"""
struct DefaultThermodynamicsBackend end
Base.broadcastable(backend::DefaultThermodynamicsBackend) = tuple(backend)

"""Condensed phase of water."""
abstract type AbstractPhase end
"""Liquid water; saturation over a plane liquid surface."""
struct Liquid <: AbstractPhase end
"""Ice; saturation over a plane ice surface."""
struct Ice <: AbstractPhase end

"""Dry-air gas constant [J kg⁻¹ K⁻¹]: `R_d(thermo, FT)`."""
function R_d end
"""Water-vapor gas constant [J kg⁻¹ K⁻¹]: `R_v(thermo, FT)`."""
function R_v end
"""Isobaric heat capacity of dry air [J kg⁻¹ K⁻¹]: `cp_d(thermo, FT)`."""
function cp_d end
"""Isobaric heat capacity of water vapor [J kg⁻¹ K⁻¹]: `cp_v(thermo, FT)`."""
function cp_v end
"""Heat capacity of liquid water [J kg⁻¹ K⁻¹]: `cp_l(thermo, FT)`."""
function cp_l end
"""Heat capacity of ice [J kg⁻¹ K⁻¹]: `cp_i(thermo, FT)`."""
function cp_i end
"""Gravitational acceleration [m s⁻²]: `grav(thermo, FT)`."""
function grav end
"""Triple-point temperature of water [K]: `T_triple(thermo, FT)`."""
function T_triple end

"""
    saturation_vapor_pressure(thermo, T, phase)

Saturation vapor pressure over a plane surface of `phase` ([`Liquid`](@ref) or [`Ice`](@ref))
[Pa] at temperature `T` [K]. A backend must satisfy `d ln e/dT = L(T) / (R_v T²)` with its own
[`latent_heat`](@ref), and give equal liquid and ice values at [`T_triple`](@ref).
"""
function saturation_vapor_pressure end

"""
    latent_heat(thermo, T, phase)

Latent heat of the transition from vapor to `phase` [J kg⁻¹] at temperature `T` [K]:
vaporization for [`Liquid`](@ref), sublimation for [`Ice`](@ref). A backend must satisfy
`dL/dT = cp_v − c`, with `c` its heat capacity of `phase`.
"""
function latent_heat end

"""Heat capacity of `phase` [J kg⁻¹ K⁻¹]."""
@inline heat_capacity(thermo, ::Liquid, ::Type{FT}) where {FT} = cp_l(thermo, FT)
@inline heat_capacity(thermo, ::Ice, ::Type{FT}) where {FT} = cp_i(thermo, FT)

@inline R_d(::DefaultThermodynamicsBackend, ::Type{FT}) where {FT} = FT(287.0)
@inline R_v(::DefaultThermodynamicsBackend, ::Type{FT}) where {FT} = FT(461.5)
@inline cp_d(::DefaultThermodynamicsBackend, ::Type{FT}) where {FT} = FT(1004.5)
@inline cp_v(::DefaultThermodynamicsBackend, ::Type{FT}) where {FT} = FT(1859)
@inline cp_l(::DefaultThermodynamicsBackend, ::Type{FT}) where {FT} = FT(4181)
@inline cp_i(::DefaultThermodynamicsBackend, ::Type{FT}) where {FT} = FT(2070)
@inline grav(::DefaultThermodynamicsBackend, ::Type{FT}) where {FT} = FT(9.81)
@inline T_triple(::DefaultThermodynamicsBackend, ::Type{FT}) where {FT} = FT(273.16)
@inline _T_0(::DefaultThermodynamicsBackend, ::Type{FT}) where {FT} = FT(273.16)
@inline _press_triple(::DefaultThermodynamicsBackend, ::Type{FT}) where {FT} = FT(611.657)
@inline _latent_heat_0(::DefaultThermodynamicsBackend, ::Liquid, ::Type{FT}) where {FT} = FT(2.5008e6)
@inline _latent_heat_0(::DefaultThermodynamicsBackend, ::Ice, ::Type{FT}) where {FT} = FT(2.8344e6)

@inline function latent_heat(backend::DefaultThermodynamicsBackend, T::FT, phase::AbstractPhase) where {FT}
    Δcp = cp_v(backend, FT) - heat_capacity(backend, phase, FT)
    return _latent_heat_0(backend, phase, FT) + Δcp * (T - _T_0(backend, FT))
end

@inline function saturation_vapor_pressure(
    backend::DefaultThermodynamicsBackend,
    T::FT,
    phase::AbstractPhase,
) where {FT}
    T_tr = T_triple(backend, FT)
    Δcp = cp_v(backend, FT) - heat_capacity(backend, phase, FT)
    L_0 = _latent_heat_0(backend, phase, FT)
    exponent = Δcp * log(T / T_tr) + (L_0 - Δcp * _T_0(backend, FT)) * (inv(T_tr) - inv(T))
    return _press_triple(backend, FT) * exp(exponent / R_v(backend, FT))
end

"""
    cp_m(thermo, q_t, q_l, q_i)

Isobaric heat capacity of moist air [J kg⁻¹ K⁻¹] per unit moist-air mass, from specific humidities.
"""
@inline function cp_m(thermo, q_t::FT, q_l::FT, q_i::FT) where {FT}
    return cp_d(thermo, FT) * (1 - q_t) + cp_v(thermo, FT) * (q_t - q_l - q_i) +
           cp_l(thermo, FT) * q_l + cp_i(thermo, FT) * q_i
end

"""
    gas_constant_air(thermo, q_t, q_l, q_i)

Gas constant of moist air [J kg⁻¹ K⁻¹], `R_d (1 − q_t) + R_v q_v`, from specific humidities.
"""
@inline function gas_constant_air(thermo, q_t::FT, q_l::FT, q_i::FT) where {FT}
    return R_d(thermo, FT) * (1 - q_t) + R_v(thermo, FT) * (q_t - q_l - q_i)
end

"""
    air_density(thermo, T, p, q_t, q_l, q_i)

Density of moist air [kg m⁻³] from the ideal-gas law, from specific humidities.
"""
@inline function air_density(thermo, T::FT, p::FT, q_t::FT, q_l::FT, q_i::FT) where {FT}
    return p / (gas_constant_air(thermo, q_t, q_l, q_i) * T)
end

"""
    PhaseSaturation{FT}

Saturation over one condensed phase at `(T, p)`: vapor pressure `e` [Pa], latent heat `L`
[J kg⁻¹] and `dL_dT`, and the saturation mixing ratio `r = ε e / (p − e)` [kg kg⁻¹ of dry air]
with `dr_dT`, `d2r_dT2`, `dr_dp`, and `d2r_dTdp`.
"""
struct PhaseSaturation{FT}
    e::FT
    L::FT
    dL_dT::FT
    r::FT
    dr_dT::FT
    d2r_dT2::FT
    dr_dp::FT
    d2r_dTdp::FT
end

"""
    saturation(thermo, T, p, phase) -> PhaseSaturation

[`PhaseSaturation`](@ref) over `phase` at temperature `T` [K] and pressure `p` [Pa]. The
temperature derivatives follow from the backend's Clausius–Clapeyron and Kirchhoff relations.
"""
@inline function saturation(thermo, T::FT, p::FT, phase::AbstractPhase) where {FT}
    e = saturation_vapor_pressure(thermo, T, phase)
    L = latent_heat(thermo, T, phase)
    dL_dT = cp_v(thermo, FT) - heat_capacity(thermo, phase, FT)
    Rv = R_v(thermo, FT)
    λ = L / (Rv * T^2)
    pe = p - e
    r = R_d(thermo, FT) / Rv * e / pe
    β = p / pe * λ
    dr_dT = r * β
    d2r_dT2 = dr_dT * (β + e * λ / pe + dL_dT / L - 2 / T)
    dr_dp = -r / pe
    d2r_dTdp = -r * λ * (p + e) / pe^2
    return PhaseSaturation{FT}(e, L, dL_dT, r, dr_dT, d2r_dT2, dr_dp, d2r_dTdp)
end
