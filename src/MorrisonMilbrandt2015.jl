"""
    MorrisonMilbrandt2015

Homogeneous-parcel condensation and deposition source solvers:

- [`MM2015PiecewiseLinear`](@ref) — piecewise-linear `S = δ/(τ Γ)` until a milestone.
- [`MM2015FixedT`](@ref) — Morrison & Milbrandt (2015) Appendix C at frozen T (event-driven C5).
- [`MM2015`](@ref) — Appendix C with residual T-updating: C6 mass, then NonEquilibrium `pθq` residuals.

Kernels take unpacked scalars `g, L, c_p, …`. Thermodynamics.jl is a weak dependency
(`params` plus variables; no stored thermodynamic state).
"""
module MorrisonMilbrandt2015

export AbstractMoistureBasis,
    DryAirMixingRatio,
    SpecificHumidity,
    MM2015State,
    MM2015Timescales,
    MM2015Forcing,
    MM2015Problem,
    MM2015PiecewiseLinear,
    MM2015FixedT,
    MM2015,
    DefaultThermodynamicsBackend,
    tendencies,
    validate

include("types.jl")
include("units.jl")
include("solvers.jl")
include("LambertW.jl")
include("thermodynamics.jl")
include("equations.jl")
include("depletion.jl")
include("schemes/MM2015PiecewiseLinear.jl")
include("schemes/MM2015FixedT.jl")
include("schemes/MM2015.jl")
include("interface.jl")
include("plotting.jl")

end # module
