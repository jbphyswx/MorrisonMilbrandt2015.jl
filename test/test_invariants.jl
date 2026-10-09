using Test: Test
using ClimaParams: ClimaParams
using Thermodynamics: Thermodynamics as TD
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

isdefined(@__MODULE__, :ParcelCorpus) || include(joinpath(@__DIR__, "corpus.jl"))

backends() = (("default", MM2015.DefaultThermodynamicsBackend()), ("Thermodynamics.jl", TD.Parameters.ThermodynamicsParameters(Float64)))

Test.@testset "Invariants over the corpus" begin
    PC = ParcelCorpus
    for (bname, thermo) in backends(),
        FT in (Float32, Float64),
        scheme in (MM2015.MM2015FixedT(), MM2015.MM2015PiecewiseLinear(), MM2015.MM2015{FT}()),
        basis in (MM2015.SpecificHumidity(), MM2015.DryAirMixingRatio())

        Test.@testset "$bname, $(nameof(typeof(scheme))), $(nameof(typeof(basis))), $FT" begin
            for case in PC.CORPUS
                problem, Δt = PC.corpus_problem(case, thermo; basis, FT)
                (; T, x_liq, x_ice) = problem.state
                S_l, S_i = MM2015.tendencies(scheme, problem, Δt)
                Test.@test isfinite(S_l) && isfinite(S_i)
                Test.@test x_liq + S_l * Δt ≥ -4 * eps(FT) * x_liq
                Test.@test x_ice + S_i * Δt ≥ -4 * eps(FT) * x_ice
                iszero(x_liq) && Test.@test S_l ≥ 0
                iszero(x_ice) && Test.@test S_i ≥ 0
                crosses_triple = scheme isa MM2015.MM2015 &&
                                 any(r -> r.event == MM2015.TripleCrossing, MM2015.trajectory(scheme, problem, Δt).segments)
                if T > MM2015.T_triple(thermo, FT) && !crosses_triple
                    Test.@test S_i ≤ 0
                    iszero(x_ice) && Test.@test S_i == 0
                end
            end
        end
    end
end
