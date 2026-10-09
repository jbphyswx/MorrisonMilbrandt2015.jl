using Test: Test
using ClimaParams: ClimaParams
using Thermodynamics: Thermodynamics as TD
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

isdefined(@__MODULE__, :ParcelCorpus) || include(joinpath(@__DIR__, "corpus.jl"))

"""`problem` with the state fields in `state` and the timescales in `timescales` replaced."""
function modified(problem; state = (;), timescales = (;))
    s, τ = problem.state, problem.timescales
    new_state = MM2015.MM2015State((get(state, name, getfield(s, name)) for name in fieldnames(MM2015.MM2015State))...)
    new_timescales = MM2015.MM2015Timescales((get(timescales, name, getfield(τ, name)) for name in fieldnames(MM2015.MM2015Timescales))...)
    return MM2015.MM2015Problem(problem.basis, problem.thermo, new_state, new_timescales, problem.forcing)
end

Test.@testset "Interface" begin
    PC = ParcelCorpus
    thermo = MM2015.DefaultThermodynamicsBackend()

    Test.@testset "validate accepts the corpus" begin
        for case in PC.CORPUS, basis in (MM2015.SpecificHumidity(), MM2015.DryAirMixingRatio())
            problem, Δt = PC.corpus_problem(case, thermo; basis)
            Test.@test MM2015.validate(problem, Δt) === nothing
        end
    end

    Test.@testset "validate rejects inputs outside the domain" begin
        problem, Δt = PC.corpus_problem(PC.corpus_case(:wbf), thermo)
        Test.@test_throws ArgumentError MM2015.validate(problem, -1.0)
        Test.@test_throws ArgumentError MM2015.validate(problem, NaN)
        Test.@test_throws ArgumentError MM2015.validate(modified(problem; state = (; x_liq = -1e-6)), Δt)
        Test.@test_throws ArgumentError MM2015.validate(modified(problem; state = (; x_ice = -1e-6)), Δt)
        Test.@test_throws ArgumentError MM2015.validate(modified(problem; state = (; x_tot = problem.state.x_liq)), Δt)
        Test.@test_throws ArgumentError MM2015.validate(modified(problem; state = (; T = -1.0)), Δt)
        Test.@test_throws ArgumentError MM2015.validate(modified(problem; state = (; p = 50.0)), Δt)
        Test.@test_throws ArgumentError MM2015.validate(modified(problem; state = (; x_tot = 1.0)), Δt)
        Test.@test_throws ArgumentError MM2015.validate(modified(problem; timescales = (; τ_liq = 0.0)), Δt)
        Test.@test_throws ArgumentError MM2015.validate(modified(problem; timescales = (; τ_ice = -1.0)), Δt)
    end

    Test.@testset "$(nameof(typeof(scheme))): closed-parcel rates convert between bases by q_d" for scheme in (MM2015.MM2015FixedT(), MM2015.MM2015PiecewiseLinear())
        for case in PC.CORPUS
            closed = merge(case, (; dq_vap_dt = 0.0))
            problem_q, Δt = PC.corpus_problem(closed, thermo)
            problem_r, _ = PC.corpus_problem(closed, thermo; basis = MM2015.DryAirMixingRatio())
            q_d = 1 - problem_q.state.x_tot
            S_q = MM2015.tendencies(scheme, problem_q, Δt)
            S_r = MM2015.tendencies(scheme, problem_r, Δt)
            scale = max(abs(S_q[1]), abs(S_q[2]))
            Test.@test abs(S_q[1] - q_d * S_r[1]) ≤ 1e-11 * scale
            Test.@test abs(S_q[2] - q_d * S_r[2]) ≤ 1e-11 * scale
        end
    end
end
