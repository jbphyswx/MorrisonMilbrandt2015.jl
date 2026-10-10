using Test: Test
using ClimaParams: ClimaParams
using Thermodynamics: Thermodynamics as TD
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

isdefined(@__MODULE__, :ParcelCorpus) || include(joinpath(@__DIR__, "corpus.jl"))

isdefined(@__MODULE__, :TestHelpers) || include(joinpath(@__DIR__, "test_helpers.jl"))
using .TestHelpers: allocated, backends

Test.@testset "MM2015PiecewiseLinear" begin
    PC = ParcelCorpus
    thermo = MM2015.DefaultThermodynamicsBackend()
    PL, FixedT = MM2015.MM2015PiecewiseLinear(), MM2015.MM2015FixedT()

    Test.@testset "first-order convergence to MM2015FixedT: $name" for name in (:warm_updraft, :wbf, :ice_only_supersaturated, :ice_subliming_liquid_growing)
        problem, _ = PC.corpus_problem(PC.corpus_case(name), thermo)
        gaps = map((0.4, 0.2, 0.1, 0.05)) do Δt
            S_pl = MM2015.tendencies(PL, problem, Δt)
            S_fixed = MM2015.tendencies(FixedT, problem, Δt)
            return abs(S_pl[1] - S_fixed[1]) + abs(S_pl[2] - S_fixed[2])
        end
        ratios = gaps[1:(end - 1)] ./ gaps[2:end]
        Test.@test all(r -> 1.9 < r < 2.1, ratios)
    end

    Test.@testset "equilibrium reached at τ, then held: $name" for (name, liquid, ice) in ((:warm_updraft, true, false), (:stiff, true, false), (:ice_only_supersaturated, false, true))
        problem, Δt = PC.corpus_problem(PC.corpus_case(name), thermo)
        tr = MM2015.trajectory(PL, problem, Δt)
        k = tr.context
        a = MM2015.ActivePhases(liquid, ice)
        (; τ, A_δ, A_δi) = MM2015.relaxation(k, a)
        first_segment = first(tr.segments)
        Test.@test first_segment.event == MM2015.Equilibrium
        Test.@test first_segment.duration == τ
        held = tr.segments[2]
        Test.@test (held.δ, held.δ_i) == (A_δ * τ, A_δi * τ)
        Test.@test MM2015.state_at(tr, Δt).δ == A_δ * τ
    end

    Test.@testset "one-phase step without events: PL − FixedT = (δ₀ − δ_eq) τ e^{−Δt/τ} / (τ_l Γ_l Δt): $name" for name in (:warm_updraft, :stiff, :sluggish)
        problem, Δt = PC.corpus_problem(PC.corpus_case(name), thermo)
        k = MM2015.coefficients(problem)
        relaxation = MM2015.relaxation(k, MM2015.ActivePhases(true, false))
        τ, A_c = relaxation.τ, relaxation.A_δ
        expected = τ < Δt ? (k.δ - A_c * τ) * τ * exp(-Δt / τ) / (k.τ_l * k.Γ_l * Δt) : NaN
        S_pl, S_fixed = MM2015.tendencies(PL, problem, Δt), MM2015.tendencies(FixedT, problem, Δt)
        if τ < Δt
            Test.@test isapprox(S_pl[1] - S_fixed[1], expected; rtol = 1e-6, atol = 1e-14 * abs(S_fixed[1]))
        else
            Test.@test S_pl[1] ≈ k.δ / (k.τ_l * k.Γ_l)
        end
    end

    Test.@testset "zero equilibrium supersaturation is an ordinary state" begin
        base = MM2015.coefficients(first(PC.corpus_problem(PC.corpus_case(:warm_updraft), thermo)))
        k = MM2015.Coefficients(base.T, base.p, base.x_l, 0.0, base.x_sl, base.x_si, 2e-6, 2e-6 + base.Δ, base.Γ_l, base.Γ_i, base.α, base.τ_l, base.τ_i, 0.0, false)
        (; Δx_l, Δx_i) = MM2015.evolve(PL, k, 60.0, MM2015.NoRecorder())
        Test.@test isfinite(Δx_l) && Δx_l > 0
        Test.@test Δx_i == 0
    end

    Test.@testset "trajectory and event bound over the corpus" begin
        for case in PC.CORPUS
            problem, Δt = PC.corpus_problem(case, thermo)
            tr = MM2015.trajectory(PL, problem, Δt)
            Test.@test tr.rates == MM2015.tendencies(PL, problem, Δt)
            Test.@test count(s -> s.event != MM2015.EndOfStep, tr.segments) ≤ MM2015.max_events(PL)
            Test.@test isapprox(sum(s -> s.duration, tr.segments), Δt; rtol = 1e-14)
        end
    end

    Test.@testset "Δt = 0" begin
        problem, _ = PC.corpus_problem(PC.corpus_case(:wbf), thermo)
        Test.@test MM2015.tendencies(PL, problem, 0.0) == (; liq = 0.0, ice = 0.0)
    end

    Test.@testset "$bname, $FT, $(nameof(typeof(basis))): inferred and allocation-free" for (bname, backend) in (
            backends()...,
            ("Thermodynamics.jl Float32 parameters", TD.Parameters.ThermodynamicsParameters(Float32)),
        ),
        FT in (Float32, Float64),
        basis in (MM2015.SpecificHumidity(), MM2015.DryAirMixingRatio())

        for name in (:wbf, :activation_ascent, :dual_depletion)
            problem, Δt = PC.corpus_problem(PC.corpus_case(name), backend; basis, FT)
            Test.@test Test.@inferred(MM2015.tendencies(PL, problem, Δt)) isa NamedTuple{(:liq, :ice), Tuple{FT, FT}}
            Test.@test allocated(MM2015.tendencies, PL, problem, Δt) == 0
            k = MM2015.coefficients(problem)
            Test.@test Test.@inferred(MM2015.tendencies(PL, k, Δt)) isa NamedTuple{(:liq, :ice), Tuple{FT, FT}}
            Test.@test allocated(MM2015.tendencies, PL, k, Δt) == 0
        end
    end
end
