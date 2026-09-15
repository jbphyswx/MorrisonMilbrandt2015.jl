module MorrisonMilbrandt2015ThermodynamicsExt

using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015
using Thermodynamics: Thermodynamics as TD

const APS = TD.Parameters.AbstractThermodynamicsParameters

"""NonEquilibrium `T` from `(p, θ_li)` at fixed `q_liq`, `q_ice`. Algebraic `pθ_li` invert."""
function MM2015.air_temperature_noneq_pθq(
    params::APS,
    p::FT,
    θ_li::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
) where {FT}
    return FT(TD.air_temperature(params, TD.pθ_li(), p, θ_li, q_tot, q_liq, q_ice))
end

function MM2015.liquid_ice_pottemp(
    params::APS,
    T::FT,
    p::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
) where {FT}
    return FT(TD.liquid_ice_pottemp_given_pressure(params, T, p, q_tot, q_liq, q_ice))
end

"""TD linearized `θ_li = (T − Lq/c_pm)/Π`, so `(∂θ_li/∂T)_{p,q} = 1/Π`."""
function MM2015.dθ_li_dT(
    params::APS,
    _T::FT,
    p::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
) where {FT}
    return FT(1) / FT(TD.exner_given_pressure(params, p, q_tot, q_liq, q_ice))
end

function MM2015.cp_m(params::APS, q_tot::FT, q_liq::FT, q_ice::FT) where {FT}
    return FT(TD.cp_m(params, q_tot, q_liq, q_ice))
end

function MM2015.air_density(
    params::APS,
    T::FT,
    p::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
) where {FT}
    return FT(TD.air_density(params, T, p, q_tot, q_liq, q_ice))
end

function MM2015.gas_constant_air(params::APS, q_tot::FT, q_liq::FT, q_ice::FT) where {FT}
    return FT(TD.gas_constant_air(params, q_tot, q_liq, q_ice))
end

function MM2015.d_gas_constant_air_dt(params::APS, dqt::FT, dql::FT, dqi::FT) where {FT}
    Rd = FT(TD.Parameters.R_d(params))
    Rv = FT(TD.Parameters.R_v(params))
    return -Rd * dqt + Rv * (dqt - dql - dqi)
end

"""Derivative of paper-based `q_sl=q_d εe_sl/(p-e_sl)`."""
function MM2015.dq_vap_sat_liq_dt(
    params::APS,
    T::FT,
    p::FT,
    q_tot::FT,
    q_sl::FT,
    dTdt::FT,
    dpdt::FT,
    dqtotdt::FT,
) where {FT}
    R_d = FT(TD.Parameters.R_d(params))
    R_v = FT(TD.Parameters.R_v(params))
    e_sl = FT(TD.saturation_vapor_pressure(params, T, TD.Liquid()))
    L_l = FT(TD.latent_heat_vapor(params, T))
    q_d = one(FT) - q_tot
    r_sl = q_sl / q_d
    drdT = (R_d / R_v) * p / (p - e_sl)^2 * e_sl * L_l / (R_v * T^2)
    drdp = -r_sl / (p - e_sl)
    return -dqtotdt * r_sl + q_d * (drdT * dTdt + drdp * dpdt)
end

"""Derivative of paper-based `q_si=q_d εe_si/(p-e_si)`."""
function MM2015.dq_vap_sat_ice_dt(
    params::APS,
    T::FT,
    p::FT,
    q_tot::FT,
    q_si::FT,
    dTdt::FT,
    dpdt::FT,
    dqtotdt::FT,
) where {FT}
    R_d = FT(TD.Parameters.R_d(params))
    R_v = FT(TD.Parameters.R_v(params))
    e_si = FT(TD.saturation_vapor_pressure(params, T, TD.Ice()))
    L_i = FT(TD.latent_heat_sublim(params, T))
    q_d = one(FT) - q_tot
    r_si = q_si / q_d
    drdT = (R_d / R_v) * p / (p - e_si)^2 * e_si * L_i / (R_v * T^2)
    drdp = -r_si / (p - e_si)
    return -dqtotdt * r_si + q_d * (drdT * dTdt + drdp * dpdt)
end

"""
TD invert `T = θ Π + Lq/c_pm` with `Π = (p/p_ref)^{R_m/c_pm}`.
`T' = θ' Π + θ Π' + (Lq' c_pm − Lq c_pm')/c_pm²`.
"""
function MM2015.dT_noneq_dt(
    params::APS,
    _T::FT,
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
    cpm = FT(TD.cp_m(params, q_tot, q_liq, q_ice))
    Rm = FT(TD.gas_constant_air(params, q_tot, q_liq, q_ice))
    Π = FT(TD.exner_given_pressure(params, p, q_tot, q_liq, q_ice))
    Lq = FT(TD.humidity_weighted_latent_heat(params, q_liq, q_ice))
    dcpm =
        FT(TD.Parameters.cp_v(params) - TD.Parameters.cp_d(params)) * dqt +
        FT(TD.Parameters.cp_l(params) - TD.Parameters.cp_v(params)) * dql +
        FT(TD.Parameters.cp_i(params) - TD.Parameters.cp_v(params)) * dqi
    dRm = MM2015.d_gas_constant_air_dt(params, dqt, dql, dqi)
    dLq = FT(TD.Parameters.LH_v0(params)) * dql + FT(TD.Parameters.LH_s0(params)) * dqi
    κ = Rm / cpm
    dκ = (dRm * cpm - Rm * dcpm) / cpm^2
    p0 = FT(TD.Parameters.p_ref_theta(params))
    dlnΠ = dκ * log(p / p0) + κ / p * dpdt
    dΠ = Π * dlnΠ
    return dθdt * Π + θ_li * dΠ + (dLq * cpm - Lq * dcpm) / cpm^2
end

@inline function MM2015.thermo_scalars(params::APS, T::FT, p::FT, q_tot::FT, q_liq::FT, q_ice::FT) where {FT}
    ρ = FT(TD.air_density(params, T, p, q_tot, q_liq, q_ice))
    e_sl = FT(TD.saturation_vapor_pressure(params, T, TD.Liquid()))
    e_si = FT(TD.saturation_vapor_pressure(params, T, TD.Ice()))
    R_d = FT(TD.Parameters.R_d(params))
    R_v = FT(TD.Parameters.R_v(params))
    ε = R_d / R_v
    q_d = one(FT) - q_tot
    r_sl = ε * e_sl / (p - e_sl)
    r_si = ε * e_si / (p - e_si)
    q_sl = q_d * r_sl
    q_si = q_d * r_si
    L_l = FT(TD.latent_heat_vapor(params, T))
    L_i = FT(TD.latent_heat_sublim(params, T))
    drsl_dT = ε * p / (p - e_sl)^2 * e_sl * L_l / (R_v * T^2)
    drsi_dT = ε * p / (p - e_si)^2 * e_si * L_i / (R_v * T^2)
    return (
        g = FT(TD.Parameters.grav(params)),
        L_i,
        L_l,
        c_p = FT(TD.cp_m(params, q_tot, q_liq, q_ice)),
        T_freeze = FT(TD.Parameters.T_triple(params)),
        dqsl_dT = q_d * drsl_dT,
        dqsi_dT = q_d * drsi_dT,
        e_sl,
        e_si,
        ρ,
        q_sl,
        q_si,
        q_vap = MM2015.vapor_specific_humidity(q_tot, q_liq, q_ice),
    )
end

end # module
