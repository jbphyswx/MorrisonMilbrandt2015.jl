include(joinpath(@__DIR__, "fixture.jl"))

function epa_closed_form_single_regime(st)
    sc = st.sc
    P = MM2015.prepare_coefficients(
        sc.g, sc.L_i, sc.L_l, sc.c_p, sc.T_freeze, sc.dqsl_dT, sc.dqsi_dT, sc.e_sl, sc.e_si,
        st.q_sl, st.q_si, st.q_liq, st.q_ice, st.δ, st.δi, st.dqvdt, st.dTdt,
    )
    A_nwL = MM2015.A_c_func_no_WBF(P.q_sl, P.g, st.w, P.c_p, P.e_sl, P.dqsl_dT, st.dqvdt, st.dTdt, st.p, st.ρ)
    A_nwI = MM2015.A_c_func_no_WBF(P.q_si, P.g, st.w, P.c_p, P.e_si, P.dqsi_dT, st.dqvdt, st.dTdt, st.p, st.ρ)
    δ = st.q_vap - P.q_sl
    δi = st.q_vap - P.q_si
    if st.q_ice == 0 && δ > 0 && st.q_liq ≥ 0
        S_ql = MM2015.S_ql_func_indiv(A_nwL, st.τ_liq, δ, st.Δt, P.Γ_l)
        return S_ql, zero(FT), true
    end
    if st.q_liq == 0 && st.T < P.T_freeze && δi > 0 && st.q_ice ≥ 0
        S_qi = MM2015.S_qi_func_indiv(A_nwI, st.τ_ice, δi, st.Δt, P.Γ_i)
        return zero(FT), S_qi, true
    end
    return zero(FT), zero(FT), false
end

Test.@testset "EP Layer 0 — primitives" begin
    A_c, τ, δ0 = FT(1e-6), FT(10), FT(5e-4)
    Test.@testset "δ_func_EP endpoints" begin
        Test.@test MM2015.δ_func_EP(A_c, τ, δ0, zero(FT)) ≈ δ0
        δ_inf = MM2015.δ_func_EP(A_c, τ, δ0, FT(1e20))
        Test.@test isapprox(δ_inf, A_c * τ; rtol = 1e-12)
        t = FT(3)
        Test.@test MM2015.δ_func_EP(A_c, τ, δ0, t) ≈ A_c * τ + (δ0 - A_c * τ) * exp(-t / τ)
    end

    Test.@testset "instantaneous source limit (Δt→0) — ice driven by δ_i" begin
        st = mm2015_build_state(;
            freezing = :BF, regime = :super, q_liq = FT(0), q_ice = FT(1e-6),
            τ_liq = FT(1e9), τ_ice = FT(10), Δt = FT(1e-12), w = FT(0),
        )
        sc = st.sc
        P = MM2015.prepare_coefficients(
            sc.g, sc.L_i, sc.L_l, sc.c_p, sc.T_freeze, sc.dqsl_dT, sc.dqsi_dT, sc.e_sl, sc.e_si,
            st.q_sl, st.q_si, st.q_liq, st.q_ice, zero(FT), zero(FT), st.dqvdt, st.dTdt,
        )
        A_nwI = MM2015.A_c_func_no_WBF(P.q_si, P.g, st.w, P.c_p, P.e_si, P.dqsi_dT, st.dqvdt, st.dTdt, st.p, st.ρ)
        δi = st.q_vap - P.q_si
        S_qi = MM2015.S_qi_func_indiv(A_nwI, st.τ_ice, δi, st.Δt, P.Γ_i)
        Test.@test isapprox(S_qi, δi / (st.τ_ice * P.Γ_i); rtol = 1e-6)
    end

    Test.@testset "A_c WBF term sign flips across freezing" begin
        for (fr, expect_gap_sign) in ((:BF, +1), (:AF, -1))
            st = mm2015_build_state(; freezing = fr, regime = :wbf, q_liq = FT(1e-4), q_ice = FT(1e-4), τ_ice = FT(10))
            sc = st.sc
            P = MM2015.prepare_coefficients(
                sc.g, sc.L_i, sc.L_l, sc.c_p, sc.T_freeze, sc.dqsl_dT, sc.dqsi_dT, sc.e_sl, sc.e_si,
                st.q_sl, st.q_si, st.q_liq, st.q_ice, zero(FT), zero(FT), st.dqvdt, st.dTdt,
            )
            A_nw = MM2015.A_c_func_no_WBF(P.q_sl, P.g, st.w, P.c_p, P.e_sl, P.dqsl_dT, st.dqvdt, st.dTdt, st.p, st.ρ)
            A_wbf = MM2015.A_c_func(st.τ_ice, P.Γ_i, P.q_sl, P.q_si, P.g, st.w, P.c_p, P.e_sl, P.L_i, P.dqsl_dT, st.dqvdt, st.dTdt, st.p, st.ρ)
            Test.@test sign(A_wbf - A_nw) == -sign(P.q_sl - P.q_si) || isapprox(P.q_sl, P.q_si; atol = 1e-14)
            Test.@test sign(P.q_sl - P.q_si) == expect_gap_sign || isapprox(P.q_sl, P.q_si; rtol = 1e-3)
        end
    end

end

Test.@testset "EP vs closed-form C5 (single-regime)" begin
    st = mm2015_build_state(;
        freezing = :BF, regime = :wbf, q_liq = FT(0), q_ice = FT(1e-6),
        τ_liq = FT(1e9), τ_ice = FT(8), Δt = FT(0.1), w = FT(0.3),
    )
    S_ref_ql, S_ref_qi, ok = epa_closed_form_single_regime(st)
    Test.@test ok
    S_ql, S_qi = mm2015_call_epa_solver(st)
    Test.@test isapprox(S_qi, S_ref_qi; rtol = 1e-10, atol = 1e-16)
    Test.@test abs(S_ql - S_ref_ql) ≤ FT(1e-12)

    st_l = mm2015_build_state(;
        freezing = :AF, regime = :super, q_liq = FT(1e-6), q_ice = FT(0),
        τ_liq = FT(8), τ_ice = FT(1e9), Δt = FT(0.5), w = FT(0),
    )
    S_ref_ql, S_ref_qi, ok = epa_closed_form_single_regime(st_l)
    Test.@test ok
    S_ql, S_qi = mm2015_call_epa_solver(st_l)
    Test.@test isapprox(S_ql, S_ref_ql; rtol = 1e-10, atol = 1e-16)
    Test.@test abs(S_qi - S_ref_qi) ≤ FT(1e-12)
end

Test.@testset "EP invariants on sat-boundary classes" begin
    for fr in (:BF, :AF), (ql, qi) in ((FT(1e-4), FT(0)), (FT(0), FT(1e-4)), (FT(1e-4), FT(1e-4)))
        st = mm2015_build_state(; freezing = fr, regime = :sub, q_liq = ql, q_ice = qi, τ_liq = FT(1e-12), τ_ice = FT(1e-12), Δt = FT(10))
        mm2015_assert_invariants(mm2015_call_epa_solver(st)..., st)
    end
    for fr in (:BF, :AF), rg in (:δ0, :δi0, :δ_eps, :δi_eps)
        st = mm2015_build_state(; freezing = fr, regime = rg, q_liq = FT(1e-6), q_ice = FT(1e-10), τ_liq = FT(1), τ_ice = FT(1e5), Δt = FT(10))
        mm2015_assert_invariants(mm2015_call_epa_solver(st)..., st)
    end
end

Test.@testset "EP contracts" begin
    st = mm2015_build_state(; Δt = FT(0), q_liq = FT(1e-4), q_ice = FT(1e-4))
    Test.@test mm2015_call_epa_solver(st) == (0.0, 0.0)
end
