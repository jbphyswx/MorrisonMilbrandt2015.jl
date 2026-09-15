include(joinpath(@__DIR__, "fixture.jl"))

Test.@testset "MM2015PiecewiseLinear milestone identity" begin
    q_sl, q_si = FT(0.008), FT(0.006)
    τ_liq, τ_ice = FT(10), FT(20)
    δ_eq, δi_eq = MM2015.get_δ_eq_point(q_sl, q_si, τ_liq, τ_ice; dδdt_no_S = FT(0))
    Test.@test δ_eq < 0
    Test.@test δi_eq > 0
    Test.@test isapprox(δi_eq - δ_eq, q_sl - q_si; rtol = 1e-12)

    regime = MM2015.Supersaturated(FT(1e-4), FT(1e-4), true)
    min_t, milestone, S_ql, S_qi, _, _ = MM2015.calculate_next_standard_milestone_time(
        regime, q_sl, q_si, FT(1e-4), FT(1e-4), FT(5e-4), FT(5e-4) + (q_sl - q_si), true, τ_liq, τ_ice,
    )
    Test.@test S_ql > 0 && S_qi > 0
    Test.@test min_t > 0
    Test.@test milestone isa MM2015.MilestoneType
end

Test.@testset "MM2015PiecewiseLinear linear S on a single-phase AF liquid parcel" begin
    st = mm2015_build_state(;
        freezing = :AF, regime = :super, q_liq = FT(1e-4), q_ice = FT(0),
        τ_liq = FT(8), τ_ice = FT(1e9), Δt = FT(0.05), w = FT(0),
    )
    S_ql, S_qi = mm2015_call_sources(MM2015.MM2015PiecewiseLinear(), st)
    Test.@test S_qi == 0 || abs(S_qi) ≤ eps(FT)
    Test.@test S_ql > 0
    mm2015_assert_invariants(S_ql, S_qi, st)
end

Test.@testset "MM2015PiecewiseLinear Δt=0" begin
    st = mm2015_build_state(; Δt = FT(0), q_liq = FT(1e-4), q_ice = FT(1e-4))
    Test.@test mm2015_call_sources(MM2015.MM2015PiecewiseLinear(), st) == (0.0, 0.0)
end
