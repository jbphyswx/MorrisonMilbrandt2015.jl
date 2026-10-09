abstract type AbstractMM2015Scheme end

"""
    Thresholds(; x_min = 0, δ_min = 0, δ_i_min = 0)

Magnitudes below which a scheme rounds the state to its limit, in the moisture basis of the problem
[kg kg⁻¹]: liquid or ice below `x_min` counts as absent and returns to the vapor during the step, and a
supersaturation over liquid below `δ_min` or over ice below `δ_i_min` in magnitude counts as
saturation. The defaults round nothing.
"""
struct Thresholds{FT}
    x_min::FT
    δ_min::FT
    δ_i_min::FT
end

Thresholds(; x_min = 0.0, δ_min = 0.0, δ_i_min = 0.0) = Thresholds(promote(x_min, δ_min, δ_i_min)...)

"""
    MM2015PiecewiseLinear(; thresholds = Thresholds())

Forward-Euler step of Appendix C from each event to the next, with rates frozen at the start of
every segment.
"""
struct MM2015PiecewiseLinear{T} <: AbstractMM2015Scheme
    thresholds::Thresholds{T}
end

MM2015PiecewiseLinear(; thresholds::Thresholds = Thresholds()) = MM2015PiecewiseLinear(thresholds)

"""
    MM2015FixedT(; thresholds = Thresholds())

Appendix C of Morrison & Milbrandt (2015): the parcel model with every coefficient frozen at the
start of the step, solved exactly.
"""
struct MM2015FixedT{T} <: AbstractMM2015Scheme
    thresholds::Thresholds{T}
end

MM2015FixedT(; thresholds::Thresholds = Thresholds()) = MM2015FixedT(thresholds)

const FrozenCoefficientScheme = Union{MM2015FixedT, MM2015PiecewiseLinear}

"""
    MM2015{FT}(; rtol, atol_q, atol_T, max_steps, root_finder, thresholds)
    MM2015(; rtol, atol_q, atol_T, max_steps, root_finder, thresholds)

The parcel model solved by an exponential Rosenbrock integrator to relative tolerance `rtol`
(default `√eps(FT)`), absolute tolerance `atol_q` [kg kg⁻¹] on liquid and ice (default
`10⁻⁴ rtol`), and absolute tolerance `atol_T` [K] on temperature (default `100 rtol`), in at most
`max_steps` steps, with events bracketed by `root_finder` and the state rounded by [`Thresholds`](@ref)
`thresholds`. `MM2015()` is `MM2015{Float64}()`.
"""
struct MM2015{FT, R <: AbstractRootFinder} <: AbstractMM2015Scheme
    rtol::FT
    atol_q::FT
    atol_T::FT
    max_steps::Int
    root_finder::R
    thresholds::Thresholds{FT}
end

function MM2015{FT}(;
    rtol = sqrt(eps(FT)),
    atol_q = FT(1e-4) * FT(rtol),
    atol_T = FT(100) * FT(rtol),
    max_steps::Int = 10_000,
    root_finder::AbstractRootFinder = BrentRootFinder(),
    thresholds::Thresholds = Thresholds(),
) where {FT}
    th = Thresholds{FT}(thresholds.x_min, thresholds.δ_min, thresholds.δ_i_min)
    return MM2015{FT, typeof(root_finder)}(FT(rtol), FT(atol_q), FT(atol_T), max_steps, root_finder, th)
end

MM2015(; kwargs...) = MM2015{Float64}(; kwargs...)

"""Mass normalization of the moisture variables of an [`MM2015Problem`](@ref)."""
abstract type AbstractMoistureBasis end

"""Water mass per unit dry-air mass."""
struct DryAirMixingRatio <: AbstractMoistureBasis end

"""Water mass per unit moist-air mass."""
struct SpecificHumidity <: AbstractMoistureBasis end

"""
    MM2015State(T, p, x_tot, x_liq, x_ice)

Parcel temperature [K], pressure [Pa], and total water, liquid, and ice in the moisture basis of
the problem.
"""
struct MM2015State{FT}
    T::FT
    p::FT
    x_tot::FT
    x_liq::FT
    x_ice::FT
end

MM2015State(T, p, x_tot, x_liq, x_ice) = MM2015State(promote(T, p, x_tot, x_liq, x_ice)...)

"""
    MM2015Timescales(τ_liq, τ_ice)

Supersaturation relaxation times of liquid and of ice [s].
"""
struct MM2015Timescales{FT}
    τ_liq::FT
    τ_ice::FT
end

MM2015Timescales(τ_liq, τ_ice) = MM2015Timescales(promote(τ_liq, τ_ice)...)

"""
    MM2015Forcing(dpdt, dTdt_external, dx_vap_dt)

External forcing, constant over the step: pressure tendency [Pa s⁻¹], temperature tendency from
radiation and mixing [K s⁻¹], and vapor tendency from mixing [s⁻¹] in the moisture basis of the
problem. A hydrostatic parcel rising at speed `w` has `dpdt = -ρ g w`.
"""
struct MM2015Forcing{FT}
    dpdt::FT
    dTdt_external::FT
    dx_vap_dt::FT
end

MM2015Forcing(dpdt, dTdt_external, dx_vap_dt) = MM2015Forcing(promote(dpdt, dTdt_external, dx_vap_dt)...)

"""
    MM2015Problem(basis, thermo, state, timescales, forcing)

Homogeneous parcel problem: moisture basis, thermodynamics backend, state, relaxation times, and
external forcing.
"""
struct MM2015Problem{FT, Basis <: AbstractMoistureBasis, Thermo}
    basis::Basis
    thermo::Thermo
    state::MM2015State{FT}
    timescales::MM2015Timescales{FT}
    forcing::MM2015Forcing{FT}
end
