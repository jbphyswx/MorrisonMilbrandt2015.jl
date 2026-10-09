"""
    tendencies(scheme, problem, Δt) -> (S_liq, S_ice)
    tendencies(scheme, k::Coefficients, Δt) -> (S_liq, S_ice)

Mean phase-change rates of liquid and ice over the step `Δt` [s], in the moisture basis of
`problem` [kg kg⁻¹ s⁻¹]. The frozen-coefficient schemes also take [`Coefficients`](@ref), for example
built from the host model's thermodynamics by [`coefficients`](@ref). `Δt = 0` gives `(0, 0)`.
"""
function tendencies(scheme::FrozenCoefficientScheme, k::Coefficients{FT}, Δt::FT) where {FT}
    iszero(Δt) && return (zero(FT), zero(FT))
    Δx_l, Δx_i = evolve(scheme, k, Δt, NoRecorder())
    return (Δx_l / Δt, Δx_i / Δt)
end

tendencies(scheme::FrozenCoefficientScheme, problem::MM2015Problem{FT}, Δt::FT) where {FT} =
    tendencies(scheme, coefficients(problem), Δt)

function tendencies(scheme::MM2015, problem::MM2015Problem{FT}, Δt::FT) where {FT}
    iszero(Δt) && return (zero(FT), zero(FT))
    Δx_l, Δx_i = evolve(scheme, problem, Δt, NoRecorder())
    return (Δx_l / Δt, Δx_i / Δt)
end

"""
    trajectory(scheme, problem, Δt)
    trajectory(scheme, k::Coefficients, Δt)

The step of [`tendencies`](@ref) with its segments recorded in a [`Trajectory`](@ref): the segments
between events of a frozen-coefficient scheme, the integrator steps of [`MM2015`](@ref).
"""
function trajectory(scheme::FrozenCoefficientScheme, k::Coefficients{FT}, Δt::FT) where {FT}
    recorder = SegmentRecorder{FT}()
    iszero(Δt) && return Trajectory(scheme, k, Δt, recorder.segments, (zero(FT), zero(FT)))
    Δx_l, Δx_i = evolve(scheme, k, Δt, recorder)
    return Trajectory(scheme, k, Δt, recorder.segments, (Δx_l / Δt, Δx_i / Δt))
end

trajectory(scheme::FrozenCoefficientScheme, problem::MM2015Problem{FT}, Δt::FT) where {FT} =
    trajectory(scheme, coefficients(problem), Δt)

function trajectory(scheme::MM2015, problem::MM2015Problem{FT}, Δt::FT) where {FT}
    recorder = SegmentRecorder{FT}()
    iszero(Δt) && return Trajectory(scheme, problem, Δt, recorder.segments, (zero(FT), zero(FT)))
    Δx_l, Δx_i = evolve(scheme, problem, Δt, recorder)
    return Trajectory(scheme, problem, Δt, recorder.segments, (Δx_l / Δt, Δx_i / Δt))
end

"""
    validate(problem, Δt)

Throw an `ArgumentError` when `problem` or the step `Δt` lies outside the domain of the schemes.
[`tendencies`](@ref) does not call it.
"""
function validate(problem::MM2015Problem, Δt)
    (; basis, thermo, state, timescales, forcing) = problem
    (; T, p, x_tot, x_liq, x_ice) = state
    values = (T, p, x_tot, x_liq, x_ice, timescales.τ_liq, timescales.τ_ice, forcing.dpdt, forcing.dTdt_external, forcing.dx_vap_dt, Δt)
    all(isfinite, values) || throw(ArgumentError("state, timescales, forcing, and Δt must be finite"))
    Δt ≥ 0 || throw(ArgumentError("Δt = $Δt must be nonnegative"))
    T > 0 || throw(ArgumentError("temperature $T K must be positive"))
    p > 0 || throw(ArgumentError("pressure $p Pa must be positive"))
    timescales.τ_liq > 0 || throw(ArgumentError("τ_liq = $(timescales.τ_liq) s must be positive"))
    timescales.τ_ice > 0 || throw(ArgumentError("τ_ice = $(timescales.τ_ice) s must be positive"))
    x_liq ≥ 0 || throw(ArgumentError("liquid $x_liq must be nonnegative"))
    x_ice ≥ 0 || throw(ArgumentError("ice $x_ice must be nonnegative"))
    x_tot ≥ x_liq + x_ice || throw(ArgumentError("vapor x_tot − x_liq − x_ice = $(x_tot - x_liq - x_ice) must be nonnegative"))
    basis isa SpecificHumidity && !(x_tot < 1) && throw(ArgumentError("total specific humidity $x_tot must be below 1"))
    for phase in (Liquid(), Ice())
        e = saturation_vapor_pressure(thermo, T, phase)
        e < p || throw(ArgumentError("saturation vapor pressure over $(nameof(typeof(phase))), $e Pa, must be below the pressure $p Pa"))
    end
    return nothing
end
