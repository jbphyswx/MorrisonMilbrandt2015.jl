using Test: Test
using ClimaParams: ClimaParams
using Thermodynamics: Thermodynamics as TD
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

isdefined(@__MODULE__, :ParcelCorpus) || include(joinpath(@__DIR__, "corpus.jl"))

allocated(f, args...) = (f(args...); @allocated f(args...))

function backends(::Type{FT}) where {FT}
    default_and_td = (("default", MM2015.DefaultThermodynamicsBackend()), ("Thermodynamics.jl Float64 parameters", TD.Parameters.ThermodynamicsParameters(Float64)))
    return FT === Float64 ? default_and_td : (default_and_td..., ("Thermodynamics.jl $FT parameters", TD.Parameters.ThermodynamicsParameters(FT)))
end

Test.@testset "Inference and allocations over the corpus" begin
    PC = ParcelCorpus
    for FT in (Float32, Float64), (bname, thermo) in backends(FT), basis in (MM2015.SpecificHumidity(), MM2015.DryAirMixingRatio())
        Test.@testset "$FT, $bname, $(nameof(typeof(basis)))" begin
            schemes = (MM2015.MM2015PiecewiseLinear(), MM2015.MM2015FixedT(), MM2015.MM2015{FT}())
            for case in PC.CORPUS
                problem, Δt = PC.corpus_problem(case, thermo; basis, FT)
                Test.@test Test.@inferred(MM2015.coefficients(problem)) isa MM2015.Coefficients{FT}
                Test.@test allocated(MM2015.coefficients, problem) == 0
                Test.@test Test.@inferred(MM2015.validate(problem, Δt)) === nothing
                Test.@test allocated(MM2015.validate, problem, Δt) == 0
                for scheme in schemes
                    Test.@test Test.@inferred(MM2015.tendencies(scheme, problem, Δt)) isa NTuple{2, FT}
                    Test.@test allocated(MM2015.tendencies, scheme, problem, Δt) == 0
                    Test.@test Test.@inferred(MM2015.trajectory(scheme, problem, Δt)).rates == MM2015.tendencies(scheme, problem, Δt)
                end
            end
            T, p, q_t, q_l, q_i = FT(261), FT(8e4), FT(4e-3), FT(2e-4), FT(1e-4)
            for phase in (MM2015.Liquid(), MM2015.Ice())
                Test.@test Test.@inferred(MM2015.saturation(thermo, T, p, phase)) isa MM2015.PhaseSaturation{FT}
                Test.@test allocated(MM2015.saturation, thermo, T, p, phase) == 0
                Test.@test Test.@inferred(MM2015.saturation_vapor_pressure(thermo, T, phase)) isa FT
                Test.@test allocated(MM2015.saturation_vapor_pressure, thermo, T, phase) == 0
                Test.@test Test.@inferred(MM2015.latent_heat(thermo, T, phase)) isa FT
                Test.@test allocated(MM2015.latent_heat, thermo, T, phase) == 0
            end
            Test.@test Test.@inferred(MM2015.cp_m(thermo, q_t, q_l, q_i)) isa FT
            Test.@test allocated(MM2015.cp_m, thermo, q_t, q_l, q_i) == 0
            Test.@test Test.@inferred(MM2015.gas_constant_air(thermo, q_t, q_l, q_i)) isa FT
            Test.@test allocated(MM2015.gas_constant_air, thermo, q_t, q_l, q_i) == 0
            Test.@test Test.@inferred(MM2015.air_density(thermo, T, p, q_t, q_l, q_i)) isa FT
            Test.@test allocated(MM2015.air_density, thermo, T, p, q_t, q_l, q_i) == 0
        end
    end
end
