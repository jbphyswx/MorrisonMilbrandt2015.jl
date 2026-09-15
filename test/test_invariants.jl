include(joinpath(@__DIR__, "fixture.jl"))

Test.@testset "moisture-source invariants (designed grid)" begin
    states_epa = mm2015_designed_states()
    states_mild = mm2015_designed_states(; include_extreme_τ = false)
    states_std = filter(st -> abs(st.T - MM2015_T_FREEZE) > FT(0.5), states_mild)

    Test.@testset "epa" begin
        n_ok = 0
        for st in states_epa
            mm2015_assert_invariants(mm2015_call_sources(MM2015.MM2015FixedT(), st)..., st)
            n_ok += 1
        end
        Test.@test n_ok == length(states_epa)
    end

    Test.@testset "standard (mild subset)" begin
        n_ok = 0
        for st in states_std
            mm2015_assert_invariants(mm2015_call_sources(MM2015.MM2015PiecewiseLinear(), st)..., st)
            n_ok += 1
        end
        Test.@test n_ok == length(states_std)
    end

    Test.@testset "mm2015 (mild subset)" begin
        n_ok = 0
        for st in states_std[1:min(end, 80)]
            mm2015_assert_invariants(mm2015_call_sources(MM2015.MM2015(), st)..., st)
            n_ok += 1
        end
        Test.@test n_ok == min(length(states_std), 80)
    end
end

Test.@testset "moisture-source Δt=0 contract" begin
    st = mm2015_build_state(; Δt = FT(0), q_liq = FT(1e-4), q_ice = FT(1e-4))
    for scheme in (MM2015.MM2015FixedT(), MM2015.MM2015(), MM2015.MM2015PiecewiseLinear())
        S_ql, S_qi = mm2015_call_sources(scheme, st)
        Test.@test S_ql == 0 && S_qi == 0
    end
end
