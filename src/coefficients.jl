"""Dry-air mass fraction `m_d / m` from the total water `x_t` of the basis."""
@inline dry_fraction(::SpecificHumidity, x_t) = 1 - x_t
@inline dry_fraction(::DryAirMixingRatio, x_t) = inv(1 + x_t)

"""Specific humidity of a water amount `x` given in the basis, with dry fraction `q_d`."""
@inline specific_humidity(::SpecificHumidity, x, _) = x
@inline specific_humidity(::DryAirMixingRatio, x, q_d) = q_d * x

"""Heat capacity per unit basis mass from the moist-air value `c_pm`."""
@inline basis_heat_capacity(::SpecificHumidity, c_pm, _) = c_pm
@inline basis_heat_capacity(::DryAirMixingRatio, c_pm, q_d) = c_pm / q_d

"""Density of the basis mass from the moist-air density `ρ`."""
@inline basis_density(::SpecificHumidity, ρ, _) = ρ
@inline basis_density(::DryAirMixingRatio, ρ, q_d) = q_d * ρ

"""
    basis_saturation(basis, s::PhaseSaturation, q_d) -> (x_s, dx_s_dT, dx_s_dp, dx_s_dxt)

Saturation humidity in the basis and its derivatives with respect to temperature, pressure, and
total water.
"""
@inline basis_saturation(::SpecificHumidity, s::PhaseSaturation, q_d) = (q_d * s.r, q_d * s.dr_dT, q_d * s.dr_dp, -s.r)
@inline basis_saturation(::DryAirMixingRatio, s::PhaseSaturation, _) = (s.r, s.dr_dT, s.dr_dp, zero(s.r))

"""
    basis_saturation_second(basis, s::PhaseSaturation, q_d) -> (x_s_TT, x_s_Tp, x_s_Txt)

Second derivatives of the saturation humidity in the basis: `∂²/∂T²`, `∂²/∂T∂p`, and `∂²/∂T∂x_t`.
"""
@inline basis_saturation_second(::SpecificHumidity, s::PhaseSaturation, q_d) = (q_d * s.d2r_dT2, q_d * s.d2r_dTdp, -s.dr_dT)
@inline basis_saturation_second(::DryAirMixingRatio, s::PhaseSaturation, _) = (s.d2r_dT2, s.d2r_dTdp, zero(s.r))

"""Derivative of the saturation humidity `x_s` in the basis with respect to the basis total water."""
@inline saturation_total_water_derivative(::SpecificHumidity, x_s, q_d) = -x_s / q_d
@inline saturation_total_water_derivative(::DryAirMixingRatio, x_s, _) = zero(x_s)

"""
    total_water_derivatives(basis, q_d, q_l, q_i) -> (dq_t, dq_l, dq_i)

Derivatives of the specific humidities with respect to the total water of the basis, at fixed liquid
and ice of the basis.
"""
@inline total_water_derivatives(::SpecificHumidity, q_d, _, _) = (one(q_d), zero(q_d), zero(q_d))
@inline total_water_derivatives(::DryAirMixingRatio, q_d, q_l, q_i) = (q_d^2, -q_d * q_l, -q_d * q_i)

"""Derivative of the basis heat capacity with respect to the basis total water, from that of `c_pm`."""
@inline basis_heat_capacity_derivative(::SpecificHumidity, dc_pm, _, _) = dc_pm
@inline basis_heat_capacity_derivative(::DryAirMixingRatio, dc_pm, c_pm, q_d) = dc_pm / q_d + c_pm

"""
    Coefficients(T, p, x_l, x_i, x_sl, x_si, δ, δ_i, Γ_l, Γ_i, α, τ_l, τ_i, A_l, below_triple)

Appendix C quantities of the parcel model at one state, in the moisture basis of the problem:
temperature `T` [K], pressure `p` [Pa], liquid `x_l` and ice `x_i`, saturation humidities `x_sl` and
`x_si`, the supersaturations `δ = x_v − x_sl` and `δ_i = x_v − x_si`, the factors `Γ_l`, `Γ_i`, `α`,
the relaxation times `τ_l` and `τ_i` [s], the external forcing `A_l` [s⁻¹] of `δ`, and whether the
parcel counts as below the triple point. The constructor adds `Δ = x_sl − x_si`, the rates
`r_l = 1/(τ_l Γ_l)`, `r_i = 1/(τ_i Γ_i)`, `k_l = 1/τ_l`, `k_i = α r_i` [s⁻¹], `k_lΔ = k_l Δ`,
`k_iΔ = k_i Δ`, and the relaxation times `τ_ice_only = 1/k_i` and `τ_both = 1/(k_l + k_i)` [s].
"""
struct Coefficients{FT}
    T::FT
    p::FT
    x_l::FT
    x_i::FT
    x_sl::FT
    x_si::FT
    δ::FT
    δ_i::FT
    Δ::FT
    Γ_l::FT
    Γ_i::FT
    α::FT
    τ_l::FT
    τ_i::FT
    A_l::FT
    below_triple::Bool
    r_l::FT
    r_i::FT
    k_l::FT
    k_i::FT
    k_lΔ::FT
    k_iΔ::FT
    τ_ice_only::FT
    τ_both::FT

    function Coefficients(T::FT, p::FT, x_l::FT, x_i::FT, x_sl::FT, x_si::FT, δ::FT, δ_i::FT, Γ_l::FT, Γ_i::FT, α::FT,
            τ_l::FT, τ_i::FT, A_l::FT, below_triple::Bool) where {FT}
        Δ = x_sl - x_si
        r_l, r_i, k_l = inv(τ_l * Γ_l), inv(τ_i * Γ_i), inv(τ_l)
        k_i = α * r_i
        return new{FT}(T, p, x_l, x_i, x_sl, x_si, δ, δ_i, Δ, Γ_l, Γ_i, α, τ_l, τ_i, A_l, below_triple, r_l, r_i, k_l, k_i,
            k_l * Δ, k_i * Δ, inv(k_i), inv(k_l + k_i))
    end
end

"""
    ThermodynamicInputs(; x_sl, x_si, dx_sl_dT, dx_si_dT, e_sl, L_v, L_s, c_pm, ρ, T_triple)

Thermodynamics of the parcel at the start of a step, from the host model, in the moisture basis of the
problem: saturation humidities over liquid and ice [kg kg⁻¹] and their temperature derivatives
[kg kg⁻¹ K⁻¹], saturation vapor pressure over liquid `e_sl` [Pa], latent heats of vaporization and
sublimation [J kg⁻¹], moist-air isobaric heat capacity `c_pm` [J kg⁻¹ K⁻¹] and density `ρ` [kg m⁻³],
and the triple-point temperature [K].
"""
struct ThermodynamicInputs{FT}
    x_sl::FT
    x_si::FT
    dx_sl_dT::FT
    dx_si_dT::FT
    e_sl::FT
    L_v::FT
    L_s::FT
    c_pm::FT
    ρ::FT
    T_triple::FT
end

ThermodynamicInputs(; x_sl, x_si, dx_sl_dT, dx_si_dT, e_sl, L_v, L_s, c_pm, ρ, T_triple) =
    ThermodynamicInputs(promote(x_sl, x_si, dx_sl_dT, dx_si_dT, e_sl, L_v, L_s, c_pm, ρ, T_triple)...)

"""
    thermodynamic_inputs(problem) -> ThermodynamicInputs

[`ThermodynamicInputs`](@ref) of `problem` from its thermodynamics backend.
"""
@inline function thermodynamic_inputs(problem::MM2015Problem{FT}) where {FT}
    (; basis, thermo, state) = problem
    (; T, p, x_tot, x_liq, x_ice) = state
    q_d = dry_fraction(basis, x_tot)
    q_t = specific_humidity(basis, x_tot, q_d)
    q_l = specific_humidity(basis, x_liq, q_d)
    q_i = specific_humidity(basis, x_ice, q_d)
    liquid = saturation(thermo, T, p, Liquid())
    ice = saturation(thermo, T, p, Ice())
    x_sl, dx_sl_dT, _, _ = basis_saturation(basis, liquid, q_d)
    x_si, dx_si_dT, _, _ = basis_saturation(basis, ice, q_d)
    return ThermodynamicInputs{FT}(
        x_sl, x_si, dx_sl_dT, dx_si_dT, liquid.e, liquid.L, ice.L,
        cp_m(thermo, q_t, q_l, q_i), air_density(thermo, T, p, q_t, q_l, q_i), T_triple(thermo, FT),
    )
end

"""
    coefficients(problem) -> Coefficients
    coefficients(basis, state, timescales, forcing, inputs::ThermodynamicInputs; δ, δ_i) -> Coefficients

[`Coefficients`](@ref) at the state of `problem`, or of a state whose thermodynamics the host model
supplies in `inputs`. The supersaturations over liquid `δ` and over ice `δ_i` [kg kg⁻¹] default to
`x_v − x_sl` and `x_v − x_si` from the state; a host model that holds them passes them.
"""
@inline coefficients(problem::MM2015Problem) =
    coefficients(problem.basis, problem.state, problem.timescales, problem.forcing, thermodynamic_inputs(problem))

@inline function coefficients(
    basis::AbstractMoistureBasis,
    state::MM2015State{FT},
    timescales::MM2015Timescales{FT},
    forcing::MM2015Forcing{FT},
    inputs::ThermodynamicInputs{FT};
    δ::FT = (state.x_tot - state.x_liq - state.x_ice) - inputs.x_sl,
    δ_i::FT = (state.x_tot - state.x_liq - state.x_ice) - inputs.x_si,
) where {FT}
    (; T, p, x_tot, x_liq, x_ice) = state
    (; x_sl, x_si, dx_sl_dT, dx_si_dT, e_sl, L_v, L_s) = inputs
    q_d = dry_fraction(basis, x_tot)
    c_p = basis_heat_capacity(basis, inputs.c_pm, q_d)
    ρ = basis_density(basis, inputs.ρ, q_d)
    dx_sl_dp = -x_sl / (p - e_sl)
    dx_sl_dxt = saturation_total_water_derivative(basis, x_sl, q_d)
    (; dpdt, dTdt_external, dx_vap_dt) = forcing
    dTdt_non_latent = dTdt_external + dpdt / (ρ * c_p)
    A_l = dx_vap_dt - dx_sl_dT * dTdt_non_latent - dx_sl_dp * dpdt - dx_sl_dxt * dx_vap_dt
    below_triple = T < inputs.T_triple || (T == inputs.T_triple && dTdt_non_latent < 0)
    return Coefficients(
        T, p, x_liq, x_ice, x_sl, x_si, δ, δ_i,
        1 + L_v / c_p * dx_sl_dT, 1 + L_s / c_p * dx_si_dT, 1 + L_s / c_p * dx_sl_dT,
        timescales.τ_liq, timescales.τ_ice, A_l, below_triple,
    )
end
