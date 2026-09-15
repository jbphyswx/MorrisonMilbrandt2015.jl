include(joinpath(@__DIR__, "fixture.jl"))

function mm2015_first_segment(st; opts = MM2015.MM2015Opts{FT}(; dqvdt = st.dqvdt, dTdt = st.dTdt))
    sc = st.sc
    θ_li = MM2015.liquid_ice_pottemp(st.thermo, st.T, st.p, st.q_tot, st.q_liq, st.q_ice)
    dθdT0 = MM2015.dθ_li_dT(st.thermo, st.T, st.p, st.q_tot, st.q_liq, st.q_ice)
    P = MM2015.prepare_coefficients(
        sc.g, sc.L_i, sc.L_l, sc.c_p, sc.T_freeze, sc.dqsl_dT, sc.dqsi_dT, sc.e_sl, sc.e_si,
        st.q_sl, st.q_si, st.q_liq, st.q_ice, st.δ, st.δi, opts.dqvdt, opts.dTdt,
    )
    A_nwL = MM2015.A_c_func_no_WBF(P.q_sl, P.g, st.w, P.c_p, P.e_sl, P.dqsl_dT, opts.dqvdt, opts.dTdt, st.p, st.ρ)
    A_nwI = MM2015.A_c_func_no_WBF(P.q_si, P.g, st.w, P.c_p, P.e_si, P.dqsi_dT, opts.dqvdt, opts.dTdt, st.p, st.ρ)
    A_wbfL = MM2015.A_c_func(
        st.τ_ice, P.Γ_i, P.q_sl, P.q_si, P.g, st.w, P.c_p, P.e_sl, P.L_i, P.dqsl_dT, opts.dqvdt, opts.dTdt, st.p, st.ρ,
    )
    A_wbfI = MM2015.A_c_func(
        st.τ_liq, P.Γ_l, P.q_si, P.q_sl, P.g, st.w, P.c_p, P.e_si, P.L_l, P.dqsi_dT, opts.dqvdt, opts.dTdt, st.p, st.ρ,
    )
    τL = MM2015.τ_func_combined(st.τ_liq, st.τ_ice, P.L_i, P.c_p, P.dqsl_dT, P.Γ_i)
    τI = MM2015.τ_func_combined(st.τ_ice, st.τ_liq, P.L_l, P.c_p, P.dqsi_dT, P.Γ_l)
    BF = st.T < P.T_freeze
    la, ia = MM2015.mm2015_active_phases(
        st.δ, st.q_liq, st.q_ice, P.q_sl, P.q_si, st.τ_liq, st.τ_ice, A_nwL, A_nwI, A_wbfL, τL, BF,
    )
    A_c, τ = if la && ia
        A_wbfL, τL
    elseif la
        A_nwL, st.τ_liq
    elseif ia
        A_nwI, st.τ_ice
    else
        A_nwL, FT(Inf)
    end
    return MM2015.mm2015_next_segment(
        st.thermo, MM2015.default_mm2015_solver, st.δ, st.q_tot, st.q_liq, st.q_ice, P.q_sl, P.q_si,
        st.τ_liq, st.τ_ice, P.Γ_l, P.Γ_i, A_nwL, A_nwI, A_wbfL, A_wbfI, τL, τI,
        BF, st.Δt, θ_li, dθdT0, opts.dTdt, st.p, st.ρ, st.T, P.g, st.w, opts.dqvdt, P.T_freeze,
    ), (; A_wbfL, A_wbfI, A_nwL, A_nwI, τL, τI, P, A_c, τ, la, ia, θ_li, dθdT0)
end

"""Chain-rule `F′` vs centered `ΔF/Δt` (check of the written derivative, not a production `F′`)."""
function mm2015_check_Fprime(st, extra, t; h = FT(1e-6))
    state_at(t_) = MM2015.residual_state_at(
        st.thermo, t_, extra.θ_li, st.p, st.ρ, extra.P.g, st.w, st.q_tot, st.q_liq, st.q_ice, st.dqvdt,
        st.δ, extra.P.q_sl, extra.P.q_si, st.τ_liq, st.τ_ice, extra.P.Γ_l, extra.P.Γ_i,
        extra.A_c, extra.τ, extra.la, extra.ia, extra.dθdT0, st.dTdt,
    )
    deriv_at(t_) = MM2015.residual_deriv_at(
        st.thermo, t_, extra.θ_li, st.p, st.ρ, extra.P.g, st.w, st.q_tot, st.q_liq, st.q_ice, st.dqvdt,
        st.δ, extra.P.q_sl, extra.P.q_si, st.τ_liq, st.τ_ice, extra.P.Γ_l, extra.P.Γ_i,
        extra.A_c, extra.τ, extra.la, extra.ia, extra.dθdT0, st.dTdt,
    )
    Fp = deriv_at(t)
    sm = state_at(t)
    Test.@test sm.q_liq > 0 || sm.q_ice > 0 || abs(st.w) > 0 || abs(st.dTdt) > 0
    fd_δ = (state_at(t + h).F_δ - state_at(t - h).F_δ) / (2h)
    fd_δi = (state_at(t + h).F_δi - state_at(t - h).F_δi) / (2h)
    fd_T = (state_at(t + h).F_T - state_at(t - h).F_T) / (2h)
    Test.@test isapprox(Fp.dF_δ, fd_δ; rtol = 2e-4, atol = 1e-12)
    Test.@test isapprox(Fp.dF_δi, fd_δi; rtol = 2e-4, atol = 1e-12)
    Test.@test isapprox(Fp.dF_T, fd_T; rtol = 2e-4, atol = 1e-10)
    return nothing
end

Test.@testset "noneq pθq invert is T→θ→T at fixed q, not saturation adjustment" begin
    backend = MM2015.DefaultThermodynamicsBackend()
    T, p = FT(261), FT(80000)
    q_tot, q_liq, q_ice = FT(0.004), FT(1e-4), FT(1e-4)
    θ = MM2015.liquid_ice_pottemp(backend, T, p, q_tot, q_liq, q_ice)
    T2 = MM2015.air_temperature_noneq_pθq(backend, p, θ, q_tot, q_liq, q_ice)
    Test.@test isapprox(T2, T; rtol = 1e-8, atol = 1e-6)

    θ_td = MM2015.liquid_ice_pottemp(MM2015_THERMO, T, p, q_tot, q_liq, q_ice)
    T_td = MM2015.air_temperature_noneq_pθq(MM2015_THERMO, p, θ_td, q_tot, q_liq, q_ice)
    Test.@test isapprox(T_td, T; rtol = 1e-8, atol = 1e-6)
end

Test.@testset "default thermo: ClimaParams table, moist cp_m, analytic ∂q*/∂T, CC at T_triple" begin
    b = MM2015.DefaultThermodynamicsBackend()
    Test.@test MM2015.grav(b) == 9.81
    Test.@test MM2015.R_d(b) == 287.0
    Test.@test MM2015.R_v(b) == 461.5
    Test.@test MM2015.cp_d(b) == 1004.5
    Test.@test MM2015.cp_v(b) == 1859
    Test.@test MM2015.cp_l(b) == 4181
    Test.@test MM2015.cp_i(b) == 2070
    Test.@test MM2015.e_ref(b) == 611.657
    Test.@test MM2015.T_triple(b) == 273.16
    Test.@test MM2015.T_freeze(b) == 273.15
    Test.@test MM2015.molmass_ratio(b) ≈ MM2015.molar_mass_water(b) / MM2015.molar_mass_dry_air(b)

    T, p = FT(261), FT(80000)
    q_tot, q_liq, q_ice = FT(0.004), FT(1e-4), FT(1e-4)
    sc = MM2015.thermo_scalars(b, T, p, q_tot, q_liq, q_ice)
    Test.@test sc.c_p ≈ MM2015.cp_m(b, q_tot, q_liq, q_ice)
    Test.@test sc.c_p != MM2015.cp_d(b)
    dT = FT(1e-4)
    fd = (
        MM2015.thermo_scalars(b, T + dT, p, q_tot, q_liq, q_ice).q_sl -
        MM2015.thermo_scalars(b, T - dT, p, q_tot, q_liq, q_ice).q_sl
    ) / (2 * dT)
    Test.@test isapprox(sc.dqsl_dT, fd; rtol = 1e-6)
end

Test.@testset "leftmost_root: soonest root, not Inf, not the later root" begin
    f = t -> (t - 1) * (t - 3)
    df = t -> 2t - 4
    tstar = MM2015.leftmost_root(f, df, 0.0, 4.0)
    Test.@test isapprox(tstar, 1.0; atol = 1e-8)
    Test.@test abs(f(tstar)) ≤ 1e-10
    Test.@test tstar < 2.0
end

Test.@testset "leftmost_root: subdivision and t_C5 both catch a fast root" begin
    # F′ = (t−1.5e-12)^2 − (0.5e-12)^2. Endpoint-only F′ probes miss both
    # turning points; deterministic subdivision and the C5 partition must not.
    c, a, C = 1.5e-12, 0.5e-12, -1e-40
    f = t -> (t - c)^3 / 3 - a^2 * (t - c) + C
    df = t -> (t - c)^2 - a^2
    t_miss = MM2015.leftmost_root(f, df, 0.0, 10.0)
    t_hit = MM2015.leftmost_root(f, df, 0.0, 10.0; t_c5 = c)
    Test.@test isfinite(t_miss)
    Test.@test isfinite(t_hit)
    Test.@test abs(f(t_miss)) ≤ 1e-20 || t_miss ≤ 3e-12
    Test.@test abs(f(t_hit)) ≤ 1e-20 || t_hit ≤ 3e-12
    Test.@test isapprox(t_miss, t_hit; atol = 3e-12)
end

Test.@testset "MM2015 residual F_δ(t⋆) ≈ 0 on noneq-updated state" begin
    st = mm2015_build_state(;
        freezing = :AF, regime = :super, q_liq = FT(1e-4), q_ice = FT(0),
        τ_liq = FT(2), τ_ice = FT(1e9), Δt = FT(20), w = FT(0),
    )
    seg, _ = mm2015_first_segment(st)
    Test.@test seg.iev == 1
    Test.@test abs(seg.F_δ) ≤ FT(1e-10)
    S_ql, S_qi = mm2015_call_sources(MM2015.MM2015(), st)
    mm2015_assert_invariants(S_ql, S_qi, st)
end

Test.@testset "MM2015 residual F_δi(t⋆) ≈ 0" begin
    st = mm2015_build_state(;
        freezing = :BF, regime = :wbf, q_liq = FT(0), q_ice = FT(1e-4),
        τ_liq = FT(1e9), τ_ice = FT(2), Δt = FT(20), w = FT(0),
    )
    seg, _ = mm2015_first_segment(st)
    Test.@test seg.iev == 2
    Test.@test abs(seg.F_δi) ≤ FT(1e-10)
end

Test.@testset "MM2015 residual F_T(t⋆) ≈ 0" begin
    # Ascent lowers p at fixed θ_li, so T falls through T_triple with no condensate events.
    st = mm2015_build_state(;
        freezing = :AF, regime = :sub, q_liq = FT(0), q_ice = FT(0),
        τ_liq = FT(1e9), τ_ice = FT(1e9), Δt = FT(80), w = FT(1),
        T = MM2015_T_FREEZE + FT(0.3),
        q_vap = FT(1e-5),
    )
    seg, _ = mm2015_first_segment(st)
    Test.@test seg.iev == 5
    Test.@test abs(seg.F_T) ≤ FT(1e-5)
end

Test.@testset "MM2015 fast-τ saturation is not jumped when Δt is large" begin
    st = mm2015_build_state(;
        freezing = :AF, regime = :super, q_liq = FT(1e-4), q_ice = FT(0),
        τ_liq = FT(1e-12), τ_ice = FT(1e9), Δt = FT(10), w = FT(0),
    )
    seg, extra = mm2015_first_segment(st)
    Test.@test isfinite(seg.dt)
    Test.@test seg.dt < FT(1e-3)
    Test.@test abs(seg.F_δ) ≤ FT(1e-8)
end

Test.@testset "analytic F′ matches centered ΔF/Δt on the residual trajectory" begin
    st_liq = mm2015_build_state(;
        freezing = :AF, regime = :super, q_liq = FT(1e-4), q_ice = FT(0),
        τ_liq = FT(2), τ_ice = FT(1e9), Δt = FT(20), w = FT(0.5), dTdt = FT(1e-3),
    )
    _, extra_liq = mm2015_first_segment(st_liq)
    mm2015_check_Fprime(st_liq, extra_liq, FT(1))

    st_wbf = mm2015_build_state(;
        freezing = :BF, regime = :wbf, q_liq = FT(1e-4), q_ice = FT(1e-4),
        τ_liq = FT(5), τ_ice = FT(8), Δt = FT(20), w = FT(0.2),
    )
    _, extra_wbf = mm2015_first_segment(st_wbf)
    mm2015_check_Fprime(st_wbf, extra_wbf, FT(1))

    backend = MM2015.DefaultThermodynamicsBackend()
    T, p, q_tot, q_liq, q_ice = FT(280), FT(80000), FT(0.004), FT(1e-4), FT(0)
    ρ = MM2015.air_density(backend, T, p, q_tot, q_liq, q_ice)
    sc = MM2015.thermo_scalars(backend, T, p, q_tot, q_liq, q_ice)
    q_sl = sc.q_sl
    q_si = sc.q_si
    q_vap = q_tot - q_liq - q_ice
    δ = q_vap - q_sl
    δi = q_vap - q_si
    st_def = (;
        thermo = backend, T, p, ρ, w = FT(0.5), τ_liq = FT(2), τ_ice = FT(1e9),
        q_tot, q_liq, q_ice, q_sl, q_si, δ, δi, dqvdt = FT(0), dTdt = FT(1e-3), sc,
        Δt = FT(20),
    )
    _, extra_def = mm2015_first_segment(st_def)
    mm2015_check_Fprime(st_def, extra_def, FT(1))
end

Test.@testset "host dTdt moves T relative to conserved-θ baseline" begin
    backend = MM2015.DefaultThermodynamicsBackend()
    T, p, ρ = FT(280), FT(80000), FT(1.0)
    q_tot, q_liq, q_ice = FT(0.004), FT(1e-4), FT(0)
    θ0 = MM2015.liquid_ice_pottemp(backend, T, p, q_tot, q_liq, q_ice)
    dθdT0 = MM2015.dθ_li_dT(backend, T, p, q_tot, q_liq, q_ice)
    t = FT(10)
    st0 = MM2015.residual_state_at(
        backend, t, θ0, p, ρ, MM2015.grav(backend), FT(0), q_tot, q_liq, q_ice, FT(0),
        FT(0), FT(0), dθdT0, FT(0),
    )
    st1 = MM2015.residual_state_at(
        backend, t, θ0, p, ρ, MM2015.grav(backend), FT(0), q_tot, q_liq, q_ice, FT(0),
        FT(0), FT(0), dθdT0, FT(1e-3),
    )
    Test.@test st1.T > st0.T
    Test.@test isapprox(st0.θ_li, θ0; rtol = 1e-12)
    Test.@test st1.θ_li ≈ θ0 + dθdT0 * FT(1e-3) * t
end

Test.@testset "default thermo zero-condensate derivative retains latent heating" begin
    backend = MM2015.DefaultThermodynamicsBackend()
    T, p, q_tot = FT(280), FT(8e4), FT(0.004)
    q_liq = q_ice = zero(FT)
    θ = MM2015.liquid_ice_pottemp(backend, T, p, q_tot, q_liq, q_ice)
    dql, dqi = FT(2e-7), zero(FT)
    analytic = MM2015.dT_noneq_dt(
        backend,
        T,
        p,
        θ,
        q_tot,
        q_liq,
        q_ice,
        zero(FT),
        zero(FT),
        zero(FT),
        dql,
        dqi,
    )
    h = FT(1e-3)
    T_h = MM2015.air_temperature_noneq_pθq(
        backend,
        p,
        θ,
        q_tot,
        dql * h,
        dqi * h,
    )
    finite_difference = (T_h - T) / h
    Test.@test isapprox(analytic, finite_difference; rtol = 2e-5, atol = 1e-9)
end

Test.@testset "thermo_scalars c_p is moist cp_m (TD backend)" begin
    T, p = FT(261), FT(80000)
    q_tot, q_liq, q_ice = FT(0.004), FT(1e-4), FT(1e-4)
    sc = MM2015.thermo_scalars(MM2015_THERMO, T, p, q_tot, q_liq, q_ice)
    Test.@test sc.c_p ≈ MM2015.cp_m(MM2015_THERMO, q_tot, q_liq, q_ice)
    Test.@test sc.c_p != TD.Parameters.cp_d(MM2015_THERMO)
end

Test.@testset "MM2015 Δt=0" begin
    st = mm2015_build_state(; Δt = FT(0), q_liq = FT(1e-4), q_ice = FT(1e-4))
    Test.@test mm2015_call_sources(MM2015.MM2015(), st) == (0.0, 0.0)
end

import RootSolvers

Test.@testset "RootSolvers extension methods are certified Phase B iterators" begin
    f = t -> (t - 1) * (t - 3)
    df = t -> 2t - 4
    for method in (
        MM2015.MM2015RSBrentSolverMethod,
        MM2015.MM2015RSRegulaFalsiSolverMethod,
        MM2015.MM2015RSBisectionSolverMethod,
        MM2015.MM2015RSNewtonADSolverMethod,
    )
        t_rs = MM2015.leftmost_root(f, df, 0.0, 4.0, MM2015.MM2015Solver{method}())
        Test.@test isapprox(t_rs, 1.0; atol = 1e-8)
        Test.@test abs(f(t_rs)) ≤ 1e-10
    end
end
