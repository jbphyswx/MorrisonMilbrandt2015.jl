using Test: Test
using ClimaParams: ClimaParams
using Thermodynamics: Thermodynamics as TD
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

isdefined(@__MODULE__, :ParcelCorpus) || include(joinpath(@__DIR__, "corpus.jl"))
isdefined(@__MODULE__, :ParcelReference) || include(joinpath(@__DIR__, "reference", "parcel_ode.jl"))

isdefined(@__MODULE__, :TestHelpers) || include(joinpath(@__DIR__, "test_helpers.jl"))
using .TestHelpers: backends

# atol_q sits above the rounding floor ε·x_v of δ = x_v − x_sl
const TIGHT = (; rtol = 1e-13, atol_q = 1e-18, atol_T = 1e-11)

"""Frozen-coefficient relaxation time and forcing of `k` for the active phases `(l, i)`."""
function frozen_τ_A(k, l, i)
    τ = inv((l ? inv(k.τ_l) : 0.0) + (i ? k.α / (k.τ_i * k.Γ_i) : 0.0))
    return τ, k.A_l - (i ? k.Δ * k.α / (k.τ_i * k.Γ_i) : 0.0)
end

"""Appendix C closed forms `(δ, x_l, x_i)` at time `t` for the active phases `(l, i)`."""
function frozen_closed_form(k, l, i, t)
    τ, A_c = frozen_τ_A(k, l, i)
    decay = -expm1(-t / τ)
    δ = A_c * τ + (k.δ - A_c * τ) * (1 - decay)
    x_l = k.x_l + (l ? (A_c * τ * t + (k.δ - A_c * τ) * τ * decay) / (k.τ_l * k.Γ_l) : 0.0)
    x_i = k.x_i + (i ? ((A_c * τ + k.Δ) * t + (k.δ - A_c * τ) * τ * decay) / (k.τ_i * k.Γ_i) : 0.0)
    return (δ, x_l, x_i)
end

"""Root of `f` in `(lo, hi)` by BigFloat bisection."""
function bigfloat_root(f, lo, hi)
    setprecision(BigFloat, 256) do
        a, b = big(lo), big(hi)
        fa = f(a)
        for _ in 1:400
            m = (a + b) / 2
            (f(m) > 0) == (fa > 0) ? (a = m) : (b = m)
        end
        return Float64((a + b) / 2)
    end
end

"""Moist enthalpy per unit basis mass relative to `T_0 = T_triple`, with the backend's latent heats at `T_0`."""
function moist_enthalpy(problem, t, u)
    (; basis, thermo, state, forcing) = problem
    FT = typeof(state.T)
    x_t = state.x_tot + forcing.dx_vap_dt * t
    x_l, x_i, T = u
    x_v = x_t - x_l - x_i
    T_0 = MM2015.T_triple(thermo, FT)
    L_v0 = MM2015.latent_heat(thermo, T_0, MM2015.Liquid())
    L_f0 = MM2015.latent_heat(thermo, T_0, MM2015.Ice()) - L_v0
    dry = basis isa MM2015.SpecificHumidity ? 1 - x_t : one(FT)
    c_p = MM2015.cp_d(thermo, FT) * dry + MM2015.cp_v(thermo, FT) * x_v + MM2015.cp_l(thermo, FT) * x_l + MM2015.cp_i(thermo, FT) * x_i
    return c_p * (T - T_0) + x_v * L_v0 - x_i * L_f0
end

Test.@testset "Reference integrator" begin
    PC = ParcelCorpus
    PR = ParcelReference
    thermo = MM2015.DefaultThermodynamicsBackend()

    Test.@testset "Radau IIA tableau: row sums and quadrature order 5" begin
        for FT in (Float64, BigFloat)
            c, A = PR.radau_tableau(FT)
            tol = 16 * eps(FT)
            for i in 1:3
                Test.@test abs(sum(A[i, :]) - c[i]) ≤ tol
            end
            for q in 1:5
                Test.@test abs(sum(A[3, j] * c[j]^(q - 1) for j in 1:3) - inv(FT(q))) ≤ tol
            end
            Test.@test abs(sum(A[3, j] * c[j]^5 for j in 1:3) - inv(FT(6))) > 1e3 * tol
        end
    end

    Test.@testset "FrozenLinear reproduces the Appendix C closed forms: $name" for (name, l, i) in ((:warm_updraft, true, false), (:sluggish, true, false), (:ice_only_supersaturated, false, true), (:band, false, true))
        problem, Δt = PC.corpus_problem(PC.corpus_case(name), thermo)
        model = PR.FrozenLinear(problem)
        sol = PR.solve(model, Δt; TIGHT...)
        Test.@test isempty(sol.events)
        expected = frozen_closed_form(model.k, l, i, Δt)
        final = sol.u[end]
        Test.@test isapprox(final[1], expected[1]; rtol = 1e-11, atol = 1e-18)
        Test.@test isapprox(final[2] - model.k.x_l, expected[2] - model.k.x_l; rtol = 1e-11, atol = 1e-20)
        Test.@test isapprox(final[3] - model.k.x_i, expected[3] - model.k.x_i; rtol = 1e-11, atol = 1e-20)
    end

    Test.@testset "event times match closed forms" begin
        problem, Δt = PC.corpus_problem(PC.corpus_case(:activation_ascent), thermo)
        model = PR.FrozenLinear(problem)
        sol = PR.solve(model, Δt; TIGHT...)
        Test.@test first(sol.events)[1] == :liquid
        Test.@test isapprox(first(sol.events)[2], -model.k.δ / model.k.A_l; rtol = 1e-12)

        for (name, l, i, index) in ((:warm_evaporation, true, false, 2), (:ice_subliming_liquid_growing, true, true, 3))
            problem, Δt = PC.corpus_problem(PC.corpus_case(name), thermo)
            model = PR.FrozenLinear(problem)
            sol = PR.solve(model, Δt; TIGHT...)
            kind, t_event = first(sol.events)
            Test.@test kind == (index == 2 ? :liquid : :ice_mass)
            k = model.k
            t_root = bigfloat_root(t -> frozen_closed_form(k, l, i, t)[index], 0, Δt)
            Test.@test isapprox(t_event, t_root; rtol = 1e-11)
        end
    end

    Test.@testset "Radau IIA converges at order 5 on the parcel model" begin
        problem, _ = PC.corpus_problem(PC.corpus_case(:warm_updraft), thermo)
        model = PR.Nonlinear(problem)
        Δt = 20.0
        reference = PR.solve(model, Δt; fixed_steps = 1024, TIGHT...).u[end]
        errors = map((8, 16, 32)) do N
            final = PR.solve(model, Δt; fixed_steps = N, TIGHT...).u[end]
            return abs(final[1] - reference[1]) + abs(final[3] - reference[3]) * 1e-6
        end
        orders = log2.(errors[1:(end - 1)] ./ errors[2:end])
        Test.@test all(o -> 4.5 < o < 5.6, orders)
    end

    Test.@testset "$bname, $(nameof(typeof(basis))): moist enthalpy is conserved in a closed parcel" for (bname, backend) in backends(), basis in (MM2015.SpecificHumidity(), MM2015.DryAirMixingRatio())
        case = merge(PC.corpus_case(:wbf), (; w = 0.0, dTdt = 0.0, dq_vap_dt = 0.0))
        problem, Δt = PC.corpus_problem(case, backend; basis)
        model = PR.Nonlinear(problem)
        sol = PR.solve(model, Δt; TIGHT...)
        Test.@test any(e -> e[1] == :liquid, sol.events)
        h0 = moist_enthalpy(problem, 0.0, sol.u[1])
        scale = MM2015.latent_heat(backend, problem.state.T, MM2015.Ice()) * (problem.state.x_liq + problem.state.x_ice)
        for (t, u) in zip(sol.t, sol.u)
            Test.@test abs(moist_enthalpy(problem, t, u) - h0) ≤ 1e-10 * scale
        end
    end

    Test.@testset "$(nameof(typeof(basis))): parcel model and frozen coefficients agree to second order in Δt: $name" for name in (:wbf, :ice_subliming_above_freezing, :warm_updraft), basis in (MM2015.SpecificHumidity(), MM2015.DryAirMixingRatio())
        problem, _ = PC.corpus_problem(PC.corpus_case(name), thermo; basis)
        gaps = map((0.4, 0.2, 0.1)) do Δt
            nonlinear = PR.solve(PR.Nonlinear(problem), Δt; TIGHT...)
            frozen = PR.solve(PR.FrozenLinear(problem), Δt; TIGHT...)
            Test.@test isempty(nonlinear.events) && isempty(frozen.events)
            S_n = PR.rates(PR.Nonlinear(problem), nonlinear, Δt)
            S_f = PR.rates(PR.FrozenLinear(problem), frozen, Δt)
            return (abs(S_n[1] - S_f[1]) + abs(S_n[2] - S_f[2])) * Δt
        end
        ratios = gaps[1:(end - 1)] ./ gaps[2:end]
        Test.@test all(r -> 3.6 < r < 4.4, ratios)
    end

    Test.@testset "BigFloat and Float64 solutions agree" begin
        problem, _ = PC.corpus_problem(PC.corpus_case(:activation_ascent), thermo)
        Δt = 5.0
        sol64 = PR.solve(PR.Nonlinear(problem), Δt; TIGHT...)
        big_problem, _ = setprecision(BigFloat, 128) do
            PC.corpus_problem(PC.corpus_case(:activation_ascent), thermo; FT = BigFloat)
        end
        sol_big = setprecision(BigFloat, 128) do
            PR.solve(PR.Nonlinear(big_problem), BigFloat(Δt); rtol = 1e-20, atol_q = 1e-28, atol_T = 1e-18)
        end
        Test.@test first(sol64.events)[1] == first(sol_big.events)[1] == :liquid
        Test.@test isapprox(first(sol64.events)[2], Float64(first(sol_big.events)[2]); rtol = 1e-11)
        Test.@test isapprox(sol64.u[end][1], Float64(sol_big.u[end][1]); rtol = 1e-10)
    end

    Test.@testset "zero-mass phase on its saturation boundary activates when its supersaturation rises" begin
        problem, _ = PC.corpus_problem(PC.corpus_case(:activation_ascent), thermo)
        base = MM2015.coefficients(problem)
        for (A_l, expected) in ((1e-6, true), (-1e-6, false))
            k = MM2015.Coefficients(base.T, base.p, 0.0, 0.0, base.x_sl, base.x_si, 0.0, base.Δ, base.Γ_l, base.Γ_i, base.α, base.τ_l, base.τ_i, A_l, false)
            active, _ = PR.activation(PR.FrozenLinear(k), 0.0, [0.0, 0.0, 0.0])
            Test.@test active.liquid == expected
            Test.@test !active.ice
        end
    end
end
