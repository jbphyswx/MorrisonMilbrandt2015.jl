using Test: Test
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

const FT = Float64

Test.@testset "C2–C7 algebra" begin
    τ_liq, τ_ice = FT(10), FT(20)
    L_i, c_p, dqsl_dT, Γ_i = FT(2.5e6), FT(1004), FT(1e-4), FT(1.25)
    τ = MM2015.τ_func(τ_liq, τ_ice, L_i, c_p, dqsl_dT, Γ_i)
    Test.@test τ > 0
    α_il = 1 + (L_i / c_p) * dqsl_dT
    Test.@test inv(τ) ≈ inv(τ_liq) + α_il / (τ_ice * Γ_i)

    A_c, δ0, Δt, Γ = FT(1e-6), FT(5e-4), FT(3), FT(1.2)
    S = MM2015.S_func_indiv(A_c, τ_liq, δ0, Δt, Γ)
    Test.@test isfinite(S)
    Test.@test MM2015.δ_func(δ0, A_c, τ_liq, zero(FT)) ≈ δ0
    Test.@test MM2015.δ_func(δ0, A_c, τ_liq, Δt) ≈
               A_c * τ_liq + (δ0 - A_c * τ_liq) * exp(-Δt / τ_liq)
end

Test.@testset "C1 reconstructed independently from C6/C7 instantaneous rates" begin
    τ_liq, τ_ice = FT(7), FT(13)
    L_v, L_s, c_p = FT(2.5e6), FT(2.834e6), FT(1010)
    dqsl_dT, dqsi_dT = FT(7e-5), FT(5e-5)
    Γ_l = 1 + (L_v / c_p) * dqsl_dT
    Γ_i = 1 + (L_s / c_p) * dqsi_dT
    α_il = 1 + (L_s / c_p) * dqsl_dT
    gap = FT(8e-4)
    δ = FT(-2e-4)
    A_ext = FT(3e-6)

    C_l = δ / (τ_liq * Γ_l)
    C_i = (δ + gap) / (τ_ice * Γ_i)
    dδ_direct = A_ext - Γ_l * C_l - α_il * C_i

    τ = inv(inv(τ_liq) + α_il / (τ_ice * Γ_i))
    A_c = A_ext - gap * α_il / (τ_ice * Γ_i)
    Test.@test dδ_direct ≈ A_c - δ / τ
end

Test.@testset "C4 ice denominator follows C1/C2/C7 reconstruction" begin
    τ_ice = FT(11)
    L_v, L_s, c_p = FT(2.5e6), FT(2.834e6), FT(1007)
    dqsl_dT, dqsi_dT = FT(9e-5), FT(4e-5)
    Γ_l = 1 + (L_v / c_p) * dqsl_dT
    Γ_i = 1 + (L_s / c_p) * dqsi_dT
    q_sl, q_si = FT(7e-3), FT(5e-3)
    g, w, e_sl, p, ρ = FT(9.81), FT(0.4), FT(700), FT(8e4), FT(1)
    dqvdt, dTdt = FT(2e-7), FT(-1e-4)

    A_ext = MM2015.A_c_func_no_WBF(
        q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ,
    )
    α_il = 1 + (L_s / c_p) * dqsl_dT
    expected = A_ext - (q_sl - q_si) * α_il / (τ_ice * Γ_i)
    implemented = MM2015.A_c_func(
        τ_ice, Γ_i, q_sl, q_si, g, w, c_p, e_sl, L_s,
        dqsl_dT, dqvdt, dTdt, p, ρ,
    )
    Test.@test implemented ≈ expected
end

function simpson_integral(f, a::FT, b::FT, n::Int = 4096)
    iseven(n) || error("Simpson interval count must be even")
    h = (b - a) / n
    total = f(a) + f(b)
    for i in 1:(n - 1)
        total += (isodd(i) ? 4 : 2) * f(a + i * h)
    end
    return total * h / 3
end

Test.@testset "C6/C7 means equal independent quadrature" begin
    τ_liq, τ_ice = FT(8), FT(17)
    Γ_l, Γ_i = FT(1.18), FT(1.31)
    α_il = FT(1.24)
    τ = inv(inv(τ_liq) + α_il / (τ_ice * Γ_i))
    A_c, δ0, gap, Δt = FT(-2e-6), FT(4e-4), FT(7e-4), FT(9)
    δ(t) = A_c * τ + (δ0 - A_c * τ) * exp(-t / τ)

    liquid_quad = simpson_integral(t -> δ(t) / (τ_liq * Γ_l), 0.0, Δt) / Δt
    ice_quad = simpson_integral(t -> (δ(t) + gap) / (τ_ice * Γ_i), 0.0, Δt) / Δt
    liquid = MM2015.S_func_no_WBF(A_c, τ, τ_liq, δ0, Δt, Γ_l)
    ice = MM2015.S_func_WBF(A_c, τ, τ_ice, δ0, Δt, Γ_i, gap, zero(FT))
    Test.@test isapprox(liquid, liquid_quad; rtol = 2e-13)
    Test.@test isapprox(ice, ice_quad; rtol = 2e-13)
end

Test.@testset "closed-parcel humidity coordinates are equivalent" begin
    q_d = FT(0.985)
    cp_q = FT(1010)
    cp_r = cp_q / q_d
    L_s = FT(2.834e6)
    dqsl_q = FT(7e-5)
    dqsi_q = FT(5e-5)
    dqsl_r = dqsl_q / q_d
    dqsi_r = dqsi_q / q_d
    τ_liq, τ_ice = FT(6), FT(14)
    Γ_i_q = 1 + (L_s / cp_q) * dqsi_q
    Γ_i_r = 1 + (L_s / cp_r) * dqsi_r
    τ_q = MM2015.τ_func(τ_liq, τ_ice, L_s, cp_q, dqsl_q, Γ_i_q)
    τ_r = MM2015.τ_func(τ_liq, τ_ice, L_s, cp_r, dqsl_r, Γ_i_r)
    Test.@test Γ_i_q ≈ Γ_i_r
    Test.@test τ_q ≈ τ_r

    A_q, δ_q, Δt, Γ_l = FT(2e-6), FT(3e-4), FT(5), FT(1.2)
    S_q = MM2015.S_func_no_WBF(A_q, τ_q, τ_liq, δ_q, Δt, Γ_l)
    S_r = MM2015.S_func_no_WBF(A_q / q_d, τ_r, τ_liq, δ_q / q_d, Δt, Γ_l)
    Test.@test S_q ≈ q_d * S_r
end

Test.@testset "units round-trip" begin
    q_tot = FT(0.01)
    q = FT(0.002)
    Test.@test MM2015.mixing_ratio_to_specific_humidity(MM2015.specific_humidity_to_mixing_ratio(q, q_tot), q_tot) ≈ q

    qdot, qtdot = FT(2e-6), FT(3e-7)
    qd = 1 - q_tot
    rdot = qdot / qd + q * qtdot / qd^2
    qdot_recovered = qd * rdot - (q / qd) * qtdot
    Test.@test qdot_recovered ≈ qdot
end
