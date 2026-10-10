module TestHelpers

using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015
using Thermodynamics: Thermodynamics as TD

"""Heap allocation in bytes after warming up the call."""
allocated(f, args...) = (f(args...); @allocated f(args...))

"""Default and Thermodynamics.jl backends used by the tests."""
backends() = (("default", MM2015.DefaultThermodynamicsBackend()), ("Thermodynamics.jl", TD.Parameters.ThermodynamicsParameters(Float64)))

"""Backends including native Float32 parameters when testing Float32."""
function backends(::Type{FT}) where {FT}
    default_and_td = (
        ("default", MM2015.DefaultThermodynamicsBackend()),
        ("Thermodynamics.jl Float64 parameters", TD.Parameters.ThermodynamicsParameters(Float64)),
    )
    return FT === Float64 ? default_and_td : (default_and_td..., ("Thermodynamics.jl $FT parameters", TD.Parameters.ThermodynamicsParameters(FT)))
end

end # module TestHelpers
