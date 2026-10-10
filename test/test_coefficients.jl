using Test: Test
using ClimaParams: ClimaParams
using Thermodynamics: Thermodynamics as TD
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

isdefined(@__MODULE__, :ParcelCorpus) || include(joinpath(@__DIR__, "corpus.jl"))
isdefined(@__MODULE__, :ParcelReference) || include(joinpath(@__DIR__, "reference", "parcel_ode.jl"))

isdefined(@__MODULE__, :TestHelpers) || include(joinpath(@__DIR__, "test_helpers.jl"))
using .TestHelpers: allocated, backends

"""Centered difference of `δ` along the parcel-model flow with active phases `a`, at `t = 0`."""
function dδdt_along_flow(model, a, h)
    u = ParcelReference.initial_state(model)
    F = ParcelReference.rhs(model, a, 0.0, u)
    δ(t, v) = ParcelReference.supersaturations(model, t, v)[1]
    return (δ(h, u .+ h .* F) - δ(-h, u .- h .* F)) / (2h)
end

"""Right side `f(δ)` of the frozen-coefficient δ equation, from `Coefficients` and flags."""
f_of(k, l, i) = k.A_l - (l ? k.δ / k.τ_l : 0.0) - (i ? k.α * k.δ_i / (k.τ_i * k.Γ_i) : 0.0)

const CASES = (:warm_updraft, :wbf, :ice_subliming_liquid_growing, :ice_only_supersaturated, :activation_moistening, :freezing_level, :float_scale_near_freeze)
const BASES = (MM2015.SpecificHumidity(), MM2015.DryAirMixingRatio())

Test.@testset "Coefficients" begin
    PC, PR = ParcelCorpus, ParcelReference

    Test.@testset "$name: δ equation from Coefficients equals dδ/dt of the parcel model" for (name, thermo) in backends()
        for case_name in CASES, basis in BASES, l in (false, true), i in (false, true)
            problem, _ = PC.corpus_problem(PC.corpus_case(case_name), thermo; basis)
            k = MM2015.coefficients(problem)
            model = PR.Nonlinear(problem)
            Test.@test isapprox(dδdt_along_flow(model, PR.Active(l, i), 1e-3), f_of(k, l, i); rtol = 1e-6, atol = 1e-15)
            δ, δ_i = PR.supersaturations(model, 0.0, PR.initial_state(model))
            Test.@test k.δ ≈ δ
            Test.@test k.δ_i ≈ δ_i
        end
    end

    Test.@testset "$name: closed parcel gives the same physics in both bases" for (name, thermo) in backends()
        for case_name in CASES
            closed = merge(PC.corpus_case(case_name), (; dq_vap_dt = 0.0))
            problem_q, _ = PC.corpus_problem(closed, thermo)
            kq = MM2015.coefficients(problem_q)
            kr = MM2015.coefficients(first(PC.corpus_problem(closed, thermo; basis = MM2015.DryAirMixingRatio())))
            q_d = 1 - problem_q.state.x_tot
            for field in (:Γ_l, :Γ_i, :α, :τ_l, :τ_i)
                Test.@test isapprox(getfield(kq, field), getfield(kr, field); rtol = 1e-13)
            end
            for field in (:δ, :δ_i, :Δ, :A_l, :x_sl, :x_si, :x_l, :x_i)
                Test.@test isapprox(getfield(kq, field), q_d * getfield(kr, field); rtol = 1e-12, atol = 1e-18)
            end
            Test.@test kq.below_triple == kr.below_triple
        end
    end

    Test.@testset "triple point: below when the non-latent temperature tendency is negative" begin
        thermo = MM2015.DefaultThermodynamicsBackend()
        T_tr = MM2015.T_triple(thermo, Float64)
        at_triple = (; name = :triple, T = T_tr, p = 8.0e4, humidity = (:RH_l, 1.0), q_l = 1e-4, q_i = 1e-4, τ_l = 10.0, τ_i = 10.0, w = 0.0, dTdt = 0.0, dq_vap_dt = 0.0, Δt = 1.0)
        for (dTdt, expected) in ((-1e-4, true), (1e-4, false))
            problem, _ = PC.corpus_problem(merge(at_triple, (; dTdt)), thermo)
            Test.@test MM2015.coefficients(problem).below_triple == expected
        end
        problem, _ = PC.corpus_problem(merge(at_triple, (; T = T_tr - 1e-3, dTdt = 1e-3)), thermo)
        Test.@test MM2015.coefficients(problem).below_triple
    end

    Test.@testset "$name, $FT, $(nameof(typeof(basis))): inferred and allocation-free" for (name, thermo) in (
            backends()...,
            ("Thermodynamics.jl Float32 parameters", TD.Parameters.ThermodynamicsParameters(Float32)),
        ),
        FT in (Float32, Float64),
        basis in BASES

        problem = MM2015.MM2015Problem(
            basis, thermo,
            MM2015.MM2015State(FT(261), FT(8e4), FT(0.002), FT(2e-4), FT(1e-4)),
            MM2015.MM2015Timescales(FT(8), FT(60)),
            MM2015.MM2015Forcing(FT(-5), FT(-2e-4), FT(1e-6)),
        )
        Test.@test Test.@inferred(MM2015.coefficients(problem)) isa MM2015.Coefficients{FT}
        Test.@test allocated(MM2015.coefficients, problem) == 0
    end
end
