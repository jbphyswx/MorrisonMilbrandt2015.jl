using Test: Test
using ClimaParams: ClimaParams
using Thermodynamics: Thermodynamics as TD
using RootSolvers: RootSolvers as RS
using LinearAlgebra: LinearAlgebra
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

isdefined(@__MODULE__, :ParcelCorpus) || include(joinpath(@__DIR__, "corpus.jl"))
isdefined(@__MODULE__, :ParcelReference) || include(joinpath(@__DIR__, "reference", "parcel_ode.jl"))

isdefined(@__MODULE__, :TestHelpers) || include(joinpath(@__DIR__, "test_helpers.jl"))
using .TestHelpers: allocated, backends

const REFERENCE_TOLERANCE = (; rtol = 1e-13, atol_q = 1e-18, atol_T = 1e-11)
const TIGHT_SCHEME = MM2015.MM2015(; rtol = 1e-10, atol_q = 1e-17, atol_T = 1e-8)
const FLOAT32_SCHEME = MM2015.MM2015{Float32}()

"""exp(M) in BigFloat by scaling and squaring of a Taylor series."""
function big_exp(M::Matrix{BigFloat})
    n = size(M, 1)
    s = max(0, ceil(Int, log2(max(LinearAlgebra.opnorm(M, 1), big(1)))) + 4)
    Z = M / big(2.0)^s
    E = Matrix{BigFloat}(LinearAlgebra.I, n, n)
    term = copy(E)
    for j in 1:400
        term = term * Z / j
        E += term
        maximum(abs, term) < big(2.0)^(-precision(BigFloat) - 8) && break
    end
    for _ in 1:s
        E = E * E
    end
    return E
end

"""(φ₀, …, φ₄)(A) from the exponential of the block matrix [[A, I, 0, 0, 0], [0, 0, I, 0, 0], …]."""
function big_phi(A::Matrix{BigFloat})
    M = zeros(BigFloat, 15, 15)
    M[1:3, 1:3] = A
    for k in 1:4
        M[(3k - 2):(3k), (3k + 1):(3k + 3)] = Matrix{BigFloat}(LinearAlgebra.I, 3, 3)
    end
    E = big_exp(M)
    return ntuple(k -> E[1:3, (3k - 2):(3k)], 5)
end

"""`A` in BigFloat with entry `m` increased by `eps()` times its magnitude."""
function perturb_entry(A::Matrix{Float64}, m::Int)
    E = big.(A)
    E[m] += big(eps()) * abs(E[m])
    return E
end

"""Mean rates of `problem` over `Δt` from the Nonlinear reference."""
function reference_rates(problem, Δt)
    model = ParcelReference.Nonlinear(problem)
    sol = ParcelReference.solve(model, Δt; REFERENCE_TOLERANCE...)
    return ParcelReference.rates(model, sol, Δt), sol
end

reference_kind(kind::MM2015.EventKind) =
    kind in (MM2015.LiquidSaturation, MM2015.LiquidExhausted) ? :liquid :
    kind == MM2015.IceExhausted ? :ice_mass : kind == MM2015.IceSaturation ? :ice_saturation : :triple

"""Moist enthalpy per unit basis mass relative to `T_triple`, with the backend's latent heats there."""
function moist_enthalpy(problem, x_l, x_i, T)
    (; basis, thermo, state) = problem
    x_v = state.x_tot - x_l - x_i
    T_0 = MM2015.T_triple(thermo, Float64)
    L_v0 = MM2015.latent_heat(thermo, T_0, MM2015.Liquid())
    L_f0 = MM2015.latent_heat(thermo, T_0, MM2015.Ice()) - L_v0
    dry = basis isa MM2015.SpecificHumidity ? 1 - state.x_tot : 1.0
    c_p = MM2015.cp_d(thermo, Float64) * dry + MM2015.cp_v(thermo, Float64) * x_v + MM2015.cp_l(thermo, Float64) * x_l + MM2015.cp_i(thermo, Float64) * x_i
    return c_p * (T - T_0) + x_v * L_v0 - x_i * L_f0
end

Test.@testset "MM2015" begin
    PC, PR = ParcelCorpus, ParcelReference
    thermo = MM2015.DefaultThermodynamicsBackend()

    Test.@testset "φ-functions against extended precision, within their conditioning" begin
        V = [1.0 0.3 -0.2; 0.1 1.0 0.4; -0.05 0.2 1.0]
        matrices = Matrix{Float64}[]
        for λ in ((-0.141, -3.4e-5, 3.3e-5), (-0.2, -2.2e-5, 0.0), (-10.0, -1e-3, 1e-4)), h in (1.0, 10.0, 100.0, 1000.0, 1e4)
            push!(matrices, h .* (V * LinearAlgebra.Diagonal(collect(λ)) / V))
        end
        push!(matrices, [-0.5 0.2 0.1; 0.0 0.0 0.0; 3.0 -1.0 -0.01] .* 20, [-0.1 2.0 0.0; -2.0 -0.1 0.0; 0.0 0.0 -5.0] .* 7, zeros(3, 3))
        for name in (:wbf, :stiff, :ice_only_supersaturated), h in (1.0, 60.0)
            problem, _ = PC.corpus_problem(PC.corpus_case(name), thermo)
            s = MM2015.parcel_state(problem, 0.0, (problem.state.x_liq, problem.state.x_ice, problem.state.T))
            (; J) = MM2015.parcel_linearization(problem, MM2015.active_phases(problem, s, MM2015.below_triple(problem, s)), s)
            push!(matrices, h .* reshape(collect(J), 3, 3))
        end
        setprecision(BigFloat, 512) do
            for A in matrices
                φ = MM2015.phi_functions(Tuple(vec(A)))
                reference = big_phi(big.(A))
                responses = map(m -> big_phi(perturb_entry(A, m)), 1:9)
                for k in 2:5
                    # componentwise condition: summed change under a relative change eps() of each entry
                    condition = maximum(sum(abs.(r[k] .- reference[k]) for r in responses))
                    Test.@test maximum(abs.(reshape(collect(φ[k]), 3, 3) .- reference[k])) ≤ 8 * (condition + eps() * maximum(abs, reference[k]))
                end
            end
        end
    end

    Test.@testset "$bname, $(nameof(typeof(basis))): linearization against finite differences and the reference right side" for (bname, backend) in backends(), basis in (MM2015.SpecificHumidity(), MM2015.DryAirMixingRatio())
        for name in (:warm_updraft, :wbf, :ice_subliming_liquid_growing, :ice_only_supersaturated, :activation_moistening, :freezing_level),
            a in (MM2015.ActivePhases(true, false), MM2015.ActivePhases(false, true), MM2015.ActivePhases(true, true))

            problem, _ = PC.corpus_problem(PC.corpus_case(name), backend; basis)
            t = 3.0
            u = (problem.state.x_liq + 1e-6, problem.state.x_ice + 1e-6, problem.state.T - 0.2)
            (; F, J, v) = MM2015.parcel_linearization(problem, a, MM2015.parcel_state(problem, t, u))
            rhs = MM2015.ParcelRHS(problem, a)
            Test.@test all(isapprox.(Tuple(F), PR.rhs(PR.Nonlinear(problem), PR.Active(a.liquid, a.ice), t, collect(u)); rtol = 1e-13, atol = 1e-25))
            for (c, step) in zip(1:3, (1e-7, 1e-7, 1e-3))
                e = ntuple(k -> k == c ? step : 0.0, 3)
                column = (Tuple(rhs(t, u .+ e)) .- Tuple(rhs(t, u .- e))) ./ (2step)
                Test.@test maximum(abs.(column .- J[(3c - 2):(3c)])) ≤ 1e-6 * maximum(abs, column)
            end
            fd_v = (Tuple(rhs(t + 1e-3, u)) .- Tuple(rhs(t - 1e-3, u))) ./ 2e-3
            Test.@test maximum(abs.(fd_v .- v)) ≤ 1e-6 * maximum(abs, fd_v)
        end
    end

    Test.@testset "exprb43 converges at order 4 with fixed steps: $name" for name in (:warm_updraft, :wbf)
        problem, _ = PC.corpus_problem(PC.corpus_case(name), thermo)
        Δt = 20.0
        u0 = (problem.state.x_liq, problem.state.x_ice, problem.state.T)
        s0 = MM2015.parcel_state(problem, 0.0, u0)
        a = MM2015.active_phases(problem, s0, MM2015.below_triple(problem, s0))
        reference = PR.solve(PR.Nonlinear(problem), Δt; REFERENCE_TOLERANCE...)
        Test.@test isempty(reference.events)
        temperature_errors = map((2, 4, 8)) do N
            h = Δt / N
            u = u0
            for n in 0:(N - 1)
                (; F, J, v) = MM2015.parcel_linearization(problem, a, MM2015.parcel_state(problem, n * h, u))
                Δu, _ = MM2015.exprb43_step(MM2015.ParcelRHS(problem, a), n * h, u, h, Tuple(F), J, v, MM2015.balancing(J))
                u = u .+ Δu
            end
            return abs(u[3] - reference.u[end][3])
        end
        orders = log2.(temperature_errors[1:(end - 1)] ./ temperature_errors[2:end])
        Test.@test all(o -> 3.8 < o < 4.4, orders)
    end

    Test.@testset "$bname, $(nameof(typeof(basis))): the parcel model solved to tolerance in Float64 and Float32" for (bname, backend) in backends(), basis in (MM2015.SpecificHumidity(), MM2015.DryAirMixingRatio())
        for case in PC.CORPUS
            problem, Δt = PC.corpus_problem(case, backend; basis)
            S_ref, reference = reference_rates(problem, Δt)
            tr = MM2015.trajectory(TIGHT_SCHEME, problem, Δt)
            scale = max(abs(S_ref[1]), abs(S_ref[2]))
            Test.@test all(abs(tr.rates[c] - S_ref[c]) ≤ 10 * TIGHT_SCHEME.rtol * scale + 10 * TIGHT_SCHEME.atol_q / Δt for c in 1:2)
            events = [(reference_kind(r.event), r.t + r.duration) for r in tr.segments if r.event != MM2015.EndOfStep]
            Test.@test length(events) == length(reference.events)
            for ((kind, t), (ref_kind, t_ref)) in zip(events, reference.events)
                Test.@test kind == ref_kind
                Test.@test isapprox(t, t_ref; rtol = 1e-7, atol = 1e-9)
            end
            final = MM2015.state_at(tr, Δt)
            s = MM2015.parcel_state(problem, Δt, (final.x_l, final.x_i, final.T))
            δ_ref, δ_i_ref = PR.supersaturations(PR.Nonlinear(problem), Δt, reference.u[end])
            (; rtol, atol_q, atol_T) = TIGHT_SCHEME
            Test.@test abs(final.δ - δ_ref) ≤ 10 * (atol_q + rtol * (s.x_v + abs(s.xsl_T) * final.T) + abs(s.xsl_T) * atol_T)
            Test.@test abs(final.δ_i - δ_i_ref) ≤ 10 * (atol_q + rtol * (s.x_v + abs(s.xsi_T) * final.T) + abs(s.xsi_T) * atol_T)
            problem32, Δt32 = PC.corpus_problem(case, backend; basis, FT = Float32)
            S32 = MM2015.tendencies(FLOAT32_SCHEME, problem32, Δt32)
            Test.@test all(abs(S32[c] - S_ref[c]) ≤ 10 * FLOAT32_SCHEME.rtol * scale + 10 * FLOAT32_SCHEME.atol_q / Δt for c in 1:2)
        end
    end

    Test.@testset "Float32 temperature increments below half an ulp accumulate" begin
        case = PC.corpus_case(:activation_cooling)
        S_ref, _ = reference_rates(PC.corpus_problem(case, thermo)...)
        problem, Δt = PC.corpus_problem(case, thermo; FT = Float32)
        tr = MM2015.trajectory(MM2015.MM2015{Float32}(; atol_q = 1.0f-13), problem, Δt)
        Test.@test any(r -> abs(case.dTdt) * r.duration < eps(problem.state.T) / 2, tr.segments)
        Test.@test abs(tr.rates[1] - S_ref[1]) ≤ 1e-3 * abs(S_ref[1])
    end

    Test.@testset "error falls with the tolerance: $name" for name in (:wbf, :activation_ascent, :ice_subliming_above_freezing)
        problem, Δt = PC.corpus_problem(PC.corpus_case(name), thermo)
        S_ref, _ = reference_rates(problem, Δt)
        errors = map((1e-6, 1e-8, 1e-10)) do rtol
            S = MM2015.tendencies(MM2015.MM2015(; rtol, atol_q = 1e-18, atol_T = 1e-10), problem, Δt)
            return max(abs(S[1] - S_ref[1]), abs(S[2] - S_ref[2])) / max(abs(S_ref[1]), abs(S_ref[2]))
        end
        Test.@test errors[1] ≤ 10 * 1e-6
        Test.@test errors[1] > errors[2] > errors[3]
    end

    Test.@testset "MM2015 and MM2015FixedT agree to second order in Δt: $name" for name in (:wbf, :ice_subliming_above_freezing, :warm_updraft)
        problem, _ = PC.corpus_problem(PC.corpus_case(name), thermo)
        gaps = map((0.4, 0.2, 0.1)) do Δt
            S = MM2015.tendencies(TIGHT_SCHEME, problem, Δt)
            S_fixed = MM2015.tendencies(MM2015.MM2015FixedT(), problem, Δt)
            return (abs(S[1] - S_fixed[1]) + abs(S[2] - S_fixed[2])) * Δt
        end
        Test.@test all(r -> 3.6 < r < 4.4, gaps[1:(end - 1)] ./ gaps[2:end])
    end

    Test.@testset "$bname, $(nameof(typeof(basis))): moist enthalpy is conserved in a closed parcel" for (bname, backend) in backends(), basis in (MM2015.SpecificHumidity(), MM2015.DryAirMixingRatio())
        case = merge(PC.corpus_case(:wbf), (; w = 0.0, dTdt = 0.0, dq_vap_dt = 0.0))
        problem, Δt = PC.corpus_problem(case, backend; basis)
        tr = MM2015.trajectory(TIGHT_SCHEME, problem, Δt)
        (; x_liq, x_ice, T) = problem.state
        h0 = moist_enthalpy(problem, x_liq, x_ice, T)
        final = MM2015.state_at(tr, Δt)
        scale = MM2015.latent_heat(backend, T, MM2015.Ice()) * (x_liq + x_ice)
        Test.@test abs(moist_enthalpy(problem, final.x_l, final.x_i, final.T) - h0) ≤ 1e-8 * scale
    end

    Test.@testset "starts exactly at the triple point decide the side from the temperature tendency" begin
        T_tr = MM2015.T_triple(thermo, Float64)
        at_triple = (; name = :triple, T = T_tr, p = 8.0e4, humidity = (:RH_l, 1.001), q_l = 1e-4, q_i = 1e-4, τ_l = 10.0, τ_i = 20.0, w = 0.0, dTdt = 0.0, dq_vap_dt = 0.0, Δt = 5.0)
        cooling, Δt = PC.corpus_problem(merge(at_triple, (; dTdt = -1e-3)), thermo)
        warming, _ = PC.corpus_problem(merge(at_triple, (; dTdt = 1e-3)), thermo)
        Test.@test MM2015.tendencies(MM2015.MM2015(), cooling, Δt)[2] > 0
        Test.@test MM2015.tendencies(MM2015.MM2015(), warming, Δt)[2] ≤ 0
    end

    Test.@testset "Brent root finder brackets sign changes, exact zeros, and runs of zeros" begin
        cases = ((s -> s - 0.3, 0.0, 1.0), (s -> s^3 - 0.001, 0.0, 2.0), (s -> s ≤ 0.25 ? 0.0 : s - 0.25, 0.0, 1.0),
            (s -> s < 0.5 ? -1.0 : 1.0, 0.0, 1.0), (s -> 0.3 < s < 0.300001 ? 0.0 : s - 0.3, 0.0, 1.0))
        for (f, lo, hi) in cases
            a, b = MM2015.bracket_root(MM2015.BrentRootFinder(), f, lo, hi, f(lo), f(hi), 1e-12)
            Test.@test f(a) ≤ 0 < f(b)
            Test.@test abs(b - a) ≤ 1e-12 + 4 * eps() * max(abs(a), abs(b))
        end
    end

    Test.@testset "the RootSolvers.jl root finder gives the same solution" begin
        scheme = MM2015.MM2015(; root_finder = MM2015.RootSolversRootFinder(RS.BrentsMethod))
        for case in PC.CORPUS
            problem, Δt = PC.corpus_problem(case, thermo)
            S_brent = MM2015.tendencies(MM2015.MM2015(), problem, Δt)
            S_rs = Test.@inferred MM2015.tendencies(scheme, problem, Δt)
            Test.@test all(isapprox(S_rs[phase], S_brent[phase]; rtol = 1e-9, atol = 1e-22) for phase in (:liq, :ice))
            Test.@test allocated(MM2015.tendencies, scheme, problem, Δt) == 0
        end
    end

    Test.@testset "limits and domain" begin
        problem, Δt = PC.corpus_problem(PC.corpus_case(:wbf), thermo)
        Test.@test_throws ErrorException MM2015.tendencies(MM2015.MM2015(; max_steps = 2, rtol = 1e-12), problem, Δt)
        problem32, Δt32 = PC.corpus_problem(PC.corpus_case(:wbf), thermo; FT = Float32)
        Test.@test_throws ArgumentError MM2015.tendencies(MM2015.MM2015(), problem32, Δt32)
        Test.@test MM2015.tendencies(MM2015.MM2015(), problem, 0.0) == (; liq = 0.0, ice = 0.0)
    end

    Test.@testset "trajectory: steps, rates, and states" begin
        for case in PC.CORPUS
            problem, Δt = PC.corpus_problem(case, thermo)
            tr = MM2015.trajectory(MM2015.MM2015(), problem, Δt)
            Test.@test tr.rates == MM2015.tendencies(MM2015.MM2015(), problem, Δt)
            Test.@test isapprox(sum(r -> r.duration, tr.segments), Δt; rtol = 1e-14)
            start = MM2015.state_at(tr, 0.0)
            Test.@test (start.x_l, start.x_i, start.T) == (problem.state.x_liq, problem.state.x_ice, problem.state.T)
            for (step, next) in zip(tr.segments[1:(end - 1)], tr.segments[2:end])
                step.event == MM2015.EndOfStep || continue
                state = MM2015.state_at(tr, prevfloat(next.t))
                Test.@test isapprox(state.T, next.T; rtol = 1e-12)
            end
        end
    end

    Test.@testset "thresholds" begin
        problem, Δt = PC.corpus_problem(PC.corpus_case(:warm_evaporation), thermo)
        small = MM2015.MM2015Problem(problem.basis, problem.thermo,
            MM2015.MM2015State(problem.state.T, problem.state.p, problem.state.x_tot - problem.state.x_liq + 1e-12, 1e-12, 0.0),
            problem.timescales, problem.forcing)
        scheme = MM2015.MM2015(; thresholds = MM2015.Thresholds(; x_min = 1e-10))
        tr = MM2015.trajectory(scheme, small, Δt)
        Test.@test tr.rates[1] == -1e-12 / Δt
        Test.@test all(s -> s.event != MM2015.LiquidExhausted, tr.segments)
        Test.@test first(tr.segments).x_l == 0
        L_v = MM2015.latent_heat(thermo, problem.state.T, MM2015.Liquid())
        c_p = MM2015.cp_m(thermo, small.state.x_tot, 1e-12, 0.0)
        Test.@test isapprox(first(tr.segments).T, problem.state.T - L_v * 1e-12 / c_p; rtol = 1e-15)
        Test.@test any(s -> s.event == MM2015.LiquidExhausted, MM2015.trajectory(MM2015.MM2015(), small, Δt).segments)

        near = PC.corpus_problem(merge(PC.corpus_case(:activation_ascent), (; humidity = (:δ, -1e-11))), thermo)[1]
        formation(tr) = first(s.t + s.duration for s in tr.segments if s.event == MM2015.LiquidSaturation)
        with = MM2015.trajectory(MM2015.MM2015(; thresholds = MM2015.Thresholds(; δ_min = 1e-10)), near, 30.0)
        without = MM2015.trajectory(MM2015.MM2015(), near, 30.0)
        Test.@test formation(with) > formation(without)
        Test.@test with.rates[1] > 0
        Test.@test MM2015.state_at(with, formation(with)).δ ≥ 1e-10
    end

    Test.@testset "$bname, $FT, $(nameof(typeof(basis))): inferred and allocation-free" for (bname, backend) in (
            backends()...,
            ("Thermodynamics.jl Float32 parameters", TD.Parameters.ThermodynamicsParameters(Float32)),
        ),
        FT in (Float32, Float64),
        basis in (MM2015.SpecificHumidity(), MM2015.DryAirMixingRatio())

        scheme = MM2015.MM2015{FT}()
        for name in (:wbf, :activation_ascent, :freezing_level)
            problem, Δt = PC.corpus_problem(PC.corpus_case(name), backend; basis, FT)
            Test.@test Test.@inferred(MM2015.tendencies(scheme, problem, Δt)) isa NamedTuple{(:liq, :ice), Tuple{FT, FT}}
            Test.@test allocated(MM2015.tendencies, scheme, problem, Δt) == 0
        end
    end
end
