include(joinpath(@__DIR__, "fixture.jl"))

function mm2015_problem(st, basis = MM2015.SpecificHumidity())
    state = if basis isa MM2015.SpecificHumidity
        MM2015.MM2015State(st.T, st.p, st.q_tot, st.q_liq, st.q_ice)
    else
        q_d = 1 - st.q_tot
        MM2015.MM2015State(
            st.T,
            st.p,
            st.q_tot / q_d,
            st.q_liq / q_d,
            st.q_ice / q_d,
        )
    end
    forcing = MM2015.MM2015Forcing(
        -st.ρ * st.sc.g * st.w,
        st.dTdt,
        st.dqvdt,
        st.dqvdt,
    )
    timescales = MM2015.MM2015Timescales(st.τ_liq, st.τ_ice)
    return MM2015.MM2015Problem(basis, st.thermo, state, timescales, forcing)
end

Test.@testset "typed problem boundary matches scalar kernels" begin
    for scheme in (MM2015.MM2015PiecewiseLinear(), MM2015.MM2015FixedT(), MM2015.MM2015())
        st = mm2015_build_state(;
            freezing = :BF,
            regime = :wbf,
            q_liq = FT(1e-4),
            q_ice = FT(1e-4),
            Δt = FT(3),
        )
        problem = mm2015_problem(st)
        actual = MM2015.tendencies(scheme, problem, st.Δt)
        Test.@test all(isfinite, actual)
        Test.@test MM2015.validate(problem, st.Δt) === nothing
    end
end

Test.@testset "closed-parcel moisture bases use one physics kernel" begin
    st = mm2015_build_state(;
        freezing = :BF,
        regime = :wbf,
        q_liq = FT(1e-4),
        q_ice = FT(1e-4),
        Δt = FT(3),
        w = FT(0),
        dqvdt = FT(0),
        dTdt = FT(0),
    )
    q_d = 1 - st.q_tot
    q_inputs = MM2015._scalar_inputs(mm2015_problem(st), st.Δt)
    r_inputs = MM2015._scalar_inputs(
        mm2015_problem(st, MM2015.DryAirMixingRatio()),
        st.Δt,
    )
    ε = TD.Parameters.R_d(st.thermo) / TD.Parameters.R_v(st.thermo)
    r_sl = ε * q_inputs.e_sl / (q_inputs.p - q_inputs.e_sl)
    Test.@test q_inputs.q_vap_eq_liq ≈ q_d * r_sl
    Test.@test r_inputs.q_vap_eq_liq ≈ r_sl
    for scheme in (MM2015.MM2015PiecewiseLinear(), MM2015.MM2015FixedT(), MM2015.MM2015())
        q_rates = MM2015.tendencies(scheme, mm2015_problem(st), st.Δt)
        r_rates = MM2015.tendencies(
            scheme,
            mm2015_problem(st, MM2015.DryAirMixingRatio()),
            st.Δt,
        )
        Test.@test q_rates[1] ≈ q_d * r_rates[1]
        Test.@test q_rates[2] ≈ q_d * r_rates[2]
    end
end

Test.@testset "validation is explicit" begin
    st = mm2015_build_state()
    problem = mm2015_problem(st)
    bad_state = MM2015.MM2015State(st.T, st.p, st.q_tot, -one(FT), st.q_ice)
    bad = MM2015.MM2015Problem(
        problem.basis,
        problem.thermo,
        bad_state,
        problem.timescales,
        problem.forcing,
    )
    Test.@test_throws ArgumentError MM2015.validate(bad, st.Δt)
    # The hot path deliberately performs no broad validation call.
    Test.@test MM2015.tendencies(MM2015.MM2015FixedT(), problem, zero(FT)) == (0, 0)
end

Test.@testset "Float32 public kernels" begin
    T = Float32
    backend = MM2015.DefaultThermodynamicsBackend()
    state = MM2015.MM2015State(T(261), T(8e4), T(0.004), T(1e-4), T(1e-4))
    timescales = MM2015.MM2015Timescales(T(8), T(12))
    forcing = MM2015.MM2015Forcing(T(0), T(0), T(0), T(0))
    problem = MM2015.MM2015Problem(
        MM2015.SpecificHumidity(),
        backend,
        state,
        timescales,
        forcing,
    )
    for scheme in (MM2015.MM2015PiecewiseLinear(), MM2015.MM2015FixedT(), MM2015.MM2015())
        rates = MM2015.tendencies(scheme, problem, T(1))
        Test.@test rates isa Tuple{T, T}
        Test.@test all(isfinite, rates)
    end
end

Test.@testset "inactive phases activate after saturation crossing" begin
    st = mm2015_build_state(;
        freezing = :AF,
        regime = :sub,
        q_liq = FT(0),
        q_ice = FT(0),
        dqvdt = FT(2e-6),
        Δt = FT(30),
    )
    problem = mm2015_problem(st)
    for scheme in (MM2015.MM2015FixedT(), MM2015.MM2015())
        rates = MM2015.tendencies(scheme, problem, st.Δt)
        Test.@test rates[1] > zero(FT)
        Test.@test rates[2] == zero(FT)
    end
end
