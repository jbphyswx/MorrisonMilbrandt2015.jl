using Test: Test
using ClimaParams: ClimaParams
using Thermodynamics: Thermodynamics as TD
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

isdefined(@__MODULE__, :ParcelCorpus) || include(joinpath(@__DIR__, "corpus.jl"))
isdefined(@__MODULE__, :ParcelReference) || include(joinpath(@__DIR__, "reference", "parcel_ode.jl"))

isdefined(@__MODULE__, :TestHelpers) || include(joinpath(@__DIR__, "test_helpers.jl"))
using .TestHelpers: allocated, backends

const REFERENCE_TOLERANCE = (; rtol = 1e-13, atol_q = 1e-18, atol_T = 1e-11)

"""Reference event symbol of a scheme event kind."""
reference_kind(kind::MM2015.EventKind) =
    kind in (MM2015.LiquidSaturation, MM2015.LiquidExhausted) ? :liquid : kind == MM2015.IceExhausted ? :ice_mass : :ice_saturation

"""Events `(kind, time)` of a recorded step, in reference symbols."""
scheme_events(trajectory) = [(reference_kind(s.event), s.t + s.duration) for s in trajectory.segments if s.event != MM2015.EndOfStep]

"""
`Coefficients` of `k` with the state fields `δ`, `x_l`, `x_i`, `A_l` replaced; `δ_i = δ + Δ` unless given.
"""
function with_state(k::MM2015.Coefficients{FT}; δ = k.δ, x_l = k.x_l, x_i = k.x_i, A_l = k.A_l, δ_i = δ + k.Δ, x_sl = k.x_sl,
        x_si = k.x_si, below_triple = k.below_triple) where {FT}
    return MM2015.Coefficients(k.T, k.p, FT(x_l), FT(x_i), FT(x_sl), FT(x_si), FT(δ), FT(δ_i), k.Γ_l, k.Γ_i, k.α, k.τ_l, k.τ_i, FT(A_l),
        below_triple)
end

"""FixedT rates of `k` over `Δt` and the rates of the FrozenLinear reference."""
function fixed_and_reference(k, Δt)
    S = MM2015.tendencies(MM2015.MM2015FixedT(), k, Δt)
    model = ParcelReference.FrozenLinear(k)
    return S, ParcelReference.rates(model, ParcelReference.solve(model, Δt; REFERENCE_TOLERANCE...), Δt)
end

rates_agree(S, S_ref, Δt) = all(abs(S[c] - S_ref[c]) ≤ 1e-10 * max(abs(S_ref[1]), abs(S_ref[2])) + 1e-17 / Δt for c in 1:2)

Test.@testset "MM2015FixedT" begin
    PC = ParcelCorpus
    PR = ParcelReference

    Test.@testset "$bname, $(nameof(typeof(basis))): exact solution of the frozen-coefficient problem" for (bname, thermo) in backends(), basis in (MM2015.SpecificHumidity(), MM2015.DryAirMixingRatio())
        for case in PC.CORPUS
            problem, Δt = PC.corpus_problem(case, thermo; basis)
            S = MM2015.tendencies(MM2015.MM2015FixedT(), problem, Δt)
            model = PR.FrozenLinear(problem)
            sol = PR.solve(model, Δt; REFERENCE_TOLERANCE...)
            Test.@test rates_agree(S, PR.rates(model, sol, Δt), Δt)
            events = scheme_events(MM2015.trajectory(MM2015.MM2015FixedT(), problem, Δt))
            Test.@test length(events) == length(sol.events)
            for ((kind, t), (ref_kind, t_ref)) in zip(events, sol.events)
                Test.@test kind == ref_kind
                Test.@test isapprox(t, t_ref; rtol = 1e-10, atol = 1e-12)
            end
        end
    end

    Test.@testset "zero-mass phases on their saturation boundaries" begin
        thermo = MM2015.DefaultThermodynamicsBackend()
        warm = MM2015.coefficients(first(PC.corpus_problem(PC.corpus_case(:activation_ascent), thermo)))
        cold = MM2015.coefficients(first(PC.corpus_problem(PC.corpus_case(:band), thermo)))
        tiny = nextfloat(0.0)
        for δ in (0.0, tiny, -tiny), A_l in (1e-6, -1e-6)
            k = with_state(warm; δ, x_l = 0.0, x_i = 0.0, A_l)
            S, S_ref = fixed_and_reference(k, 30.0)
            Test.@test rates_agree(S, S_ref, 30.0)
            Test.@test (S[1] > 0) == (A_l > 0)
        end
        for offset in (0.0, tiny, -tiny), A_l in (5e-6, -5e-6)
            k = with_state(cold; δ = -cold.Δ + offset, δ_i = offset, x_l = 0.0, x_i = 0.0, A_l)
            S, S_ref = fixed_and_reference(k, 60.0)
            Test.@test rates_agree(S, S_ref, 60.0)
            Test.@test (S[2] > 0) == (A_l > 0)
        end
    end

    Test.@testset "ice with mass above the triple point, on its saturation boundary" begin
        thermo = MM2015.DefaultThermodynamicsBackend()
        warm = MM2015.coefficients(first(PC.corpus_problem(PC.corpus_case(:ice_subliming_liquid_growing), thermo)))
        for A_l in (-3e-4, 3e-4)
            k = with_state(warm; δ = -warm.Δ, δ_i = 0.0, x_l = 1e-4, x_i = 5e-5, A_l)
            S, S_ref = fixed_and_reference(k, 20.0)
            Test.@test rates_agree(S, S_ref, 20.0)
            Test.@test S[2] ≤ 0
        end
    end

    Test.@testset "regression classes" begin
        thermo = MM2015.DefaultThermodynamicsBackend()
        tend(name) = MM2015.tendencies(MM2015.MM2015FixedT(), PC.corpus_problem(PC.corpus_case(name), thermo)...)
        S_band = tend(:band)
        Test.@test S_band[1] == 0
        Test.@test S_band[2] > 0
        Test.@test tend(:massless_ice_wbf)[2] > 0
        Test.@test tend(:ice_subliming_above_freezing)[2] < 0
        for name in (:activation_moistening, :activation_ascent, :activation_cooling)
            problem, Δt = PC.corpus_problem(PC.corpus_case(name), thermo)
            S = MM2015.tendencies(MM2015.MM2015FixedT(), problem, Δt)
            Test.@test S[1] > 0
            Test.@test S[2] == 0
            Test.@test first(MM2015.trajectory(MM2015.MM2015FixedT(), problem, Δt).segments).event == MM2015.LiquidSaturation
        end
    end

    Test.@testset "trajectory: segments, rates, and states" begin
        thermo = MM2015.DefaultThermodynamicsBackend()
        for case in PC.CORPUS
            problem, Δt = PC.corpus_problem(case, thermo)
            tr = MM2015.trajectory(MM2015.MM2015FixedT(), problem, Δt)
            Test.@test tr.rates == MM2015.tendencies(MM2015.MM2015FixedT(), problem, Δt)
            Test.@test count(s -> s.event != MM2015.EndOfStep, tr.segments) ≤ MM2015.max_events(MM2015.MM2015FixedT())
            Test.@test last(tr.segments).event == MM2015.EndOfStep
            Test.@test isapprox(sum(s -> s.duration, tr.segments), Δt; rtol = 1e-14)
            k = tr.context
            start = MM2015.state_at(tr, 0.0)
            Test.@test (start.δ, start.δ_i, start.x_l, start.x_i, start.T) == (k.δ, k.δ_i, k.x_l, k.x_i, k.T)
            final = MM2015.state_at(tr, Δt)
            Test.@test isapprox(final.x_l - k.x_l, tr.rates[1] * Δt; rtol = 1e-12, atol = 1e-22)
            Test.@test isapprox(final.x_i - k.x_i, tr.rates[2] * Δt; rtol = 1e-12, atol = 1e-22)
        end
    end

    Test.@testset "$(nameof(typeof(basis))): host thermodynamics give the coefficients of the backend" for basis in (MM2015.SpecificHumidity(), MM2015.DryAirMixingRatio())
        thermo = MM2015.DefaultThermodynamicsBackend()
        for case in PC.CORPUS
            problem, Δt = PC.corpus_problem(case, thermo; basis)
            liquid = MM2015.saturation(thermo, problem.state.T, problem.state.p, MM2015.Liquid())
            ice = MM2015.saturation(thermo, problem.state.T, problem.state.p, MM2015.Ice())
            q_d = basis isa MM2015.SpecificHumidity ? 1 - problem.state.x_tot : inv(1 + problem.state.x_tot)
            w = basis isa MM2015.SpecificHumidity ? q_d : 1.0
            q = (problem.state.x_tot, problem.state.x_liq, problem.state.x_ice) .* (basis isa MM2015.SpecificHumidity ? 1.0 : q_d)
            inputs = MM2015.ThermodynamicInputs(;
                x_sl = w * liquid.r, x_si = w * ice.r, dx_sl_dT = w * liquid.dr_dT, dx_si_dT = w * ice.dr_dT, e_sl = liquid.e,
                L_v = liquid.L, L_s = ice.L, c_pm = MM2015.cp_m(thermo, q...),
                ρ = MM2015.air_density(thermo, problem.state.T, problem.state.p, q...), T_triple = MM2015.T_triple(thermo, Float64),
            )
            host = MM2015.coefficients(basis, problem.state, problem.timescales, problem.forcing, inputs)
            Test.@test all(getfield(host, f) === getfield(MM2015.coefficients(problem), f) for f in fieldnames(MM2015.Coefficients))
            Test.@test MM2015.tendencies(MM2015.MM2015FixedT(), host, Δt) === MM2015.tendencies(MM2015.MM2015FixedT(), problem, Δt)
            given = MM2015.coefficients(basis, problem.state, problem.timescales, problem.forcing, inputs;
                δ = nextfloat(host.δ), δ_i = prevfloat(host.δ_i))
            Test.@test (given.δ, given.δ_i) == (nextfloat(host.δ), prevfloat(host.δ_i))
        end
    end

    Test.@testset "a parcel without condensate that stays subsaturated has zero rates in one segment" begin
        thermo = MM2015.DefaultThermodynamicsBackend()
        problem = MM2015.MM2015Problem(MM2015.SpecificHumidity(), thermo, MM2015.MM2015State(280.0, 9.0e4, 5.0e-3, 0.0, 0.0),
            MM2015.MM2015Timescales(5.0, 1000.0), MM2015.MM2015Forcing(-1.0, 0.0, 0.0))
        k = MM2015.coefficients(problem)
        Test.@test k.δ + k.A_l * 60.0 < 0
        Test.@test MM2015.tendencies(MM2015.MM2015FixedT(), k, 60.0) == (; liq = 0.0, ice = 0.0)
        tr = MM2015.trajectory(MM2015.MM2015FixedT(), k, 60.0)
        Test.@test length(tr.segments) == 1
        Test.@test only(tr.segments).active == MM2015.ActivePhases(false, false)
        long = 1.01 * -k.δ / k.A_l
        Test.@test MM2015.tendencies(MM2015.MM2015FixedT(), k, long)[1] > 0
    end

    Test.@testset "thresholds" begin
        thermo = MM2015.DefaultThermodynamicsBackend()
        Δt = 60.0
        evaporating = MM2015.coefficients(first(PC.corpus_problem(PC.corpus_case(:warm_evaporation), thermo)))
        k = with_state(evaporating; x_l = 1e-12)
        scheme = MM2015.MM2015FixedT(; thresholds = MM2015.Thresholds(; x_min = 1e-10))
        Test.@test MM2015.tendencies(scheme, k, Δt) == (; liq = -k.x_l / Δt, ice = 0.0)
        Test.@test all(s -> s.event != MM2015.LiquidExhausted, MM2015.trajectory(scheme, k, Δt).segments)
        Test.@test MM2015.tendencies(MM2015.MM2015FixedT(), k, Δt) == (; liq = -k.x_l / Δt, ice = 0.0)
        Test.@test any(s -> s.event == MM2015.LiquidExhausted, MM2015.trajectory(MM2015.MM2015FixedT(), k, Δt).segments)

        near = with_state(evaporating; δ = 1e-10, δ_i = 1e-10 + evaporating.Δ, x_l = 0.0)
        snapped = MM2015.trajectory(MM2015.MM2015FixedT(; thresholds = MM2015.Thresholds(; δ_min = 1e-9)), near, Δt)
        Test.@test first(snapped.segments).δ == 0
        Test.@test first(MM2015.trajectory(MM2015.MM2015FixedT(), near, Δt).segments).δ == 1e-10

        cold = MM2015.coefficients(first(PC.corpus_problem(PC.corpus_case(:band), thermo)))
        near_ice = with_state(cold; δ_i = -1e-10, δ = -1e-10 - cold.Δ)
        Test.@test first(MM2015.trajectory(MM2015.MM2015FixedT(; thresholds = MM2015.Thresholds(; δ_i_min = 1e-9)), near_ice, Δt).segments).δ_i == 0
    end

    Test.@testset "Float32: the ice rate keeps the precision of a small δ_i given with the coefficients" begin
        thermo = MM2015.DefaultThermodynamicsBackend()
        base = MM2015.coefficients(first(PC.corpus_problem(PC.corpus_case(:band), thermo; FT = Float32)))
        Δt = 60.0f0
        for δ_i in (-2.0f-9, 3.0f-9, -1.0f-8), A_l in (0.0f0, 1.0f-11)
            k32 = with_state(base; δ_i, δ = δ_i - base.Δ, x_l = 0.0f0, A_l)
            S32 = MM2015.tendencies(MM2015.MM2015FixedT(), k32, Δt)
            S_big = setprecision(BigFloat, 256) do
                kbig = MM2015.Coefficients(big(k32.T), big(k32.p), big(0.0), big(k32.x_i), big(k32.x_sl), big(k32.x_si),
                    big(δ_i) - (big(k32.x_sl) - big(k32.x_si)), big(δ_i), big(k32.Γ_l), big(k32.Γ_i), big(k32.α), big(k32.τ_l),
                    big(k32.τ_i), big(A_l), k32.below_triple)
                MM2015.tendencies(MM2015.MM2015FixedT(), kbig, big(Δt))
            end
            Test.@test S32[1] == 0
            Test.@test abs(S32[2] - S_big[2]) ≤ 64 * eps(Float32) * abs(S_big[2])
        end
    end

    Test.@testset "Δt = 0" begin
        problem, _ = PC.corpus_problem(PC.corpus_case(:wbf), MM2015.DefaultThermodynamicsBackend())
        Test.@test MM2015.tendencies(MM2015.MM2015FixedT(), problem, 0.0) == (; liq = 0.0, ice = 0.0)
    end

    Test.@testset "$bname, $FT, $(nameof(typeof(basis))): inferred and allocation-free" for (bname, thermo) in (
            backends()...,
            ("Thermodynamics.jl Float32 parameters", TD.Parameters.ThermodynamicsParameters(Float32)),
        ),
        FT in (Float32, Float64),
        basis in (MM2015.SpecificHumidity(), MM2015.DryAirMixingRatio())

        for name in (:wbf, :activation_ascent, :dual_depletion)
            problem, Δt = PC.corpus_problem(PC.corpus_case(name), thermo; basis, FT)
            Test.@test Test.@inferred(MM2015.tendencies(MM2015.MM2015FixedT(), problem, Δt)) isa NamedTuple{(:liq, :ice), Tuple{FT, FT}}
            Test.@test allocated(MM2015.tendencies, MM2015.MM2015FixedT(), problem, Δt) == 0
            k = MM2015.coefficients(problem)
            Test.@test Test.@inferred(MM2015.tendencies(MM2015.MM2015FixedT(), k, Δt)) isa NamedTuple{(:liq, :ice), Tuple{FT, FT}}
            Test.@test allocated(MM2015.tendencies, MM2015.MM2015FixedT(), k, Δt) == 0
        end
    end
end
