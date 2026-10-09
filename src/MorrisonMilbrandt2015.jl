"""
    MorrisonMilbrandt2015

Condensation, evaporation, deposition, and sublimation of a homogeneous air parcel over a time step,
after Morrison & Milbrandt (2015), Appendix C. [`tendencies`](@ref) returns the mean rates of an
[`MM2015Problem`](@ref) with one of three schemes:

- [`MM2015PiecewiseLinear`](@ref): forward-Euler segments with rates frozen at each segment start.
- [`MM2015FixedT`](@ref): the exact solution of the frozen-coefficient problem.
- [`MM2015`](@ref): the parcel model solved numerically, with temperature evolving.

Thermodynamics.jl parameter sets serve as backends through a package extension.
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
    Thresholds,
    ThermodynamicInputs,
    coefficients,
    tendencies,
    trajectory,
    validate

include("solvers.jl")
include("types.jl")
include("LambertW.jl")
include("thermodynamics.jl")
include("coefficients.jl")
include("events.jl")
include("equations.jl")
include("depletion.jl")
include("recorder.jl")
include("exprb.jl")
include("schemes/MM2015FixedT.jl")
include("schemes/MM2015PiecewiseLinear.jl")
include("schemes/MM2015.jl")
include("interface.jl")
include("plotting.jl")
include("precompile_workload.jl")

end # module
