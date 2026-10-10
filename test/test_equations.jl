using Test: Test
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

isdefined(@__MODULE__, :ParcelCorpus) || include(joinpath(@__DIR__, "corpus.jl"))

"""Composite Simpson rule on `n` intervals."""
function simpson(f, a, b, n = 4096)
    h = (b - a) / n
    return h / 3 * (f(a) + f(b) + sum((isodd(i) ? 4 : 2) * f(a + i * h) for i in 1:(n - 1)))
end

const ACTIVE_SETS = (MM2015.ActivePhases(true, false), MM2015.ActivePhases(false, true), MM2015.ActivePhases(true, true), MM2015.ActivePhases(false, false))

Test.@testset "Appendix C equations on frozen coefficients" begin
    thermo = MM2015.DefaultThermodynamicsBackend()
    for name in (:wbf, :ice_only_supersaturated, :ice_subliming_liquid_growing, :activation_ascent), a in ACTIVE_SETS
        problem, _ = ParcelCorpus.corpus_problem(ParcelCorpus.corpus_case(name), thermo)
        k = MM2015.coefficients(problem)
        (; inv_τ, τ, A_δ, A_δi) = MM2015.relaxation(k, a)
        t = 7.0
        evolution(s) = MM2015.frozen_evolution(k.δ, k.δ_i, inv_τ, τ, A_δ, A_δi, s)
        δ_of(s) = evolution(s).δ
        δ_i_of(s) = evolution(s).δ_i

        Test.@testset "$name, $a: relaxation time is the inverse of the relaxation rate" begin
            Test.@test iszero(inv_τ) ? isinf(τ) : isapprox(τ * inv_τ, 1; rtol = 4eps())
        end

        Test.@testset "$name, $a: δ_i − δ stays Δ" begin
            for s in (0.0, 0.5, t, 100.0)
                Test.@test isapprox(δ_i_of(s) - δ_of(s), k.Δ; rtol = 1e-12, atol = 1e-12 * abs(k.δ))
            end
        end

        Test.@testset "$name, $a: C5 integrals equal quadrature" begin
            (; I, I_i) = evolution(t)
            Test.@test isapprox(I, simpson(δ_of, 0.0, t); rtol = 1e-12, atol = 1e-20)
            Test.@test isapprox(I_i, simpson(δ_i_of, 0.0, t); rtol = 1e-12, atol = 1e-20)
        end

        Test.@testset "$name, $a: C1 is the time derivative of C5" begin
            h = 1e-4
            rate = MM2015.supersaturation_tendency(k, a, δ_of(t), δ_i_of(t))
            Test.@test isapprox((δ_of(t + h) - δ_of(t - h)) / (2h), rate; rtol = 1e-7, atol = 1e-16)
            Test.@test isapprox((δ_i_of(t + h) - δ_i_of(t - h)) / (2h), rate; rtol = 1e-7, atol = 1e-16)
        end

        Test.@testset "$name, $a: C6 and C7 increments integrate the rates" begin
            h = 1e-4
            increments(s) = MM2015.condensate_increments(k, a, evolution(s).I, evolution(s).I_i)
            rate_l = (increments(t + h).Δx_l - increments(t - h).Δx_l) / (2h)
            rate_i = (increments(t + h).Δx_i - increments(t - h).Δx_i) / (2h)
            Test.@test isapprox(rate_l, a.liquid ? δ_of(t) / (k.τ_l * k.Γ_l) : 0.0; rtol = 1e-7, atol = 1e-18)
            Test.@test isapprox(rate_i, a.ice ? δ_i_of(t) / (k.τ_i * k.Γ_i) : 0.0; rtol = 1e-7, atol = 1e-18)
        end

        Test.@testset "$name, $a: saturation time inverts C5" begin
            for (x_0, A, x_of) in ((k.δ, A_δ, δ_of), (k.δ_i, A_δi, δ_i_of))
                t_hit = MM2015.saturation_time(x_0, inv_τ, τ, A)
                if isfinite(t_hit)
                    Test.@test t_hit > 0
                    Test.@test abs(x_of(t_hit)) ≤ 1e-12 * max(abs(x_0), abs(k.Δ))
                else
                    Test.@test iszero(inv_τ) ? iszero(A) || -x_0 / A ≤ 0 : x_0 * A * τ ≥ 0 || iszero(x_0)
                end
            end
        end
    end
end
