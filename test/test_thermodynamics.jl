using Test: Test
using ClimaParams: ClimaParams
using Thermodynamics: Thermodynamics as TD
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

const TP = TD.Parameters
const PHASES = (MM2015.Liquid(), MM2015.Ice())
const TEMPERATURES = (200.0, 230.0, 250.0, 273.16, 285.0, 310.0)
const PRESSURES = (2.5e4, 5.0e4, 8.0e4, 1.0e5)

centered(f, x, h) = (f(x + h) - f(x - h)) / (2h)

allocated(f, args...) = (f(args...); @allocated f(args...))

function backends()
    return (
        ("default", MM2015.DefaultThermodynamicsBackend()),
        ("Thermodynamics.jl", TP.ThermodynamicsParameters(Float64)),
    )
end

Test.@testset "Thermodynamics backends" begin
    Test.@testset "default backend constants are the ClimaParams values" begin
        td = TP.ThermodynamicsParameters(Float64)
        default = MM2015.DefaultThermodynamicsBackend()
        for (ours, theirs) in (
            (MM2015.R_d, TP.R_d), (MM2015.R_v, TP.R_v), (MM2015.cp_d, TP.cp_d), (MM2015.cp_v, TP.cp_v),
            (MM2015.cp_l, TP.cp_l), (MM2015.cp_i, TP.cp_i), (MM2015.grav, TP.grav), (MM2015.T_triple, TP.T_triple),
        )
            Test.@test ours(default, Float64) == theirs(td)
        end
    end

    Test.@testset "$name: Clausius–Clapeyron with the backend's own latent heat" for (name, thermo) in backends()
        R_v = MM2015.R_v(thermo, Float64)
        for T in TEMPERATURES, phase in PHASES
            dlne_dT = centered(x -> log(MM2015.saturation_vapor_pressure(thermo, x, phase)), T, 1e-3)
            L = MM2015.latent_heat(thermo, T, phase)
            Test.@test isapprox(dlne_dT, L / (R_v * T^2); rtol = 1e-8)
        end
    end

    Test.@testset "$name: Kirchhoff with the backend's own heat capacities" for (name, thermo) in backends()
        for T in TEMPERATURES, phase in PHASES
            dL_dT = centered(x -> MM2015.latent_heat(thermo, x, phase), T, 1e-2)
            Δcp = MM2015.cp_v(thermo, Float64) - MM2015.heat_capacity(thermo, phase, Float64)
            Test.@test isapprox(dL_dT, Δcp; rtol = 1e-9)
        end
    end

    Test.@testset "$name: liquid and ice saturation agree at the triple point" for (name, thermo) in backends()
        T_tr = MM2015.T_triple(thermo, Float64)
        e_l = MM2015.saturation_vapor_pressure(thermo, T_tr, MM2015.Liquid())
        e_i = MM2015.saturation_vapor_pressure(thermo, T_tr, MM2015.Ice())
        Test.@test isapprox(e_l, e_i; rtol = 4 * eps(Float64))
    end

    Test.@testset "default backend: saturation derivatives match BigFloat differences" begin
        thermo = MM2015.DefaultThermodynamicsBackend()
        h = big(1e-15)
        for T0 in TEMPERATURES, p0 in PRESSURES, phase in PHASES
            T, p = big(T0), big(p0)
            s = MM2015.saturation(thermo, T, p, phase)
            r(T_, p_) = MM2015.saturation(thermo, T_, p_, phase).r
            dr_dT(T_, p_) = MM2015.saturation(thermo, T_, p_, phase).dr_dT
            Test.@test isapprox(s.dr_dT, centered(x -> r(x, p), T, h); rtol = 1e-20)
            Test.@test isapprox(s.d2r_dT2, centered(x -> dr_dT(x, p), T, h); rtol = 1e-20)
            Test.@test isapprox(s.dr_dp, centered(x -> r(T, x), p, h * p); rtol = 1e-20)
            Test.@test isapprox(s.d2r_dTdp, centered(x -> dr_dT(T, x), p, h * p); rtol = 1e-20)
            Test.@test isapprox(s.dL_dT, centered(x -> MM2015.latent_heat(thermo, x, phase), T, h); rtol = 1e-20)
        end
    end

    Test.@testset "Thermodynamics.jl backend: saturation derivatives match differences" begin
        thermo = TP.ThermodynamicsParameters(Float64)
        for T in TEMPERATURES, p in PRESSURES, phase in PHASES
            s = MM2015.saturation(thermo, T, p, phase)
            r(T_, p_) = MM2015.saturation(thermo, T_, p_, phase).r
            dr_dT(T_, p_) = MM2015.saturation(thermo, T_, p_, phase).dr_dT
            Test.@test isapprox(s.dr_dT, centered(x -> r(x, p), T, 1e-3); rtol = 1e-7)
            Test.@test isapprox(s.d2r_dT2, centered(x -> dr_dT(x, p), T, 1e-3); rtol = 1e-6)
            Test.@test isapprox(s.dr_dp, centered(x -> r(T, x), p, 1e-4 * p); rtol = 1e-7)
            Test.@test isapprox(s.d2r_dTdp, centered(x -> dr_dT(T, x), p, 1e-4 * p); rtol = 1e-6)
        end
    end

    Test.@testset "moist-air properties match Thermodynamics.jl" begin
        td = TP.ThermodynamicsParameters(Float64)
        for (q_t, q_l, q_i) in ((0.0, 0.0, 0.0), (0.012, 3e-4, 0.0), (0.004, 1e-4, 2e-4), (0.02, 0.0, 1e-3))
            Test.@test isapprox(MM2015.cp_m(td, q_t, q_l, q_i), TD.cp_m(td, q_t, q_l, q_i); rtol = 1e-14)
            Test.@test isapprox(MM2015.gas_constant_air(td, q_t, q_l, q_i), TD.gas_constant_air(td, q_t, q_l, q_i); rtol = 1e-14)
            Test.@test isapprox(
                MM2015.air_density(td, 265.0, 7.5e4, q_t, q_l, q_i),
                TD.air_density(td, 265.0, 7.5e4, q_t, q_l, q_i);
                rtol = 1e-14,
            )
        end
    end

    Test.@testset "$name, $FT: inferred and allocation-free" for (name, thermo) in (
            backends()...,
            ("Thermodynamics.jl Float32 parameters", TP.ThermodynamicsParameters(Float32)),
        ),
        FT in (Float32, Float64)

        T, p = FT(255), FT(7e4)
        q_t, q_l, q_i = FT(0.004), FT(1e-4), FT(5e-5)
        for phase in PHASES
            Test.@test Test.@inferred(MM2015.saturation_vapor_pressure(thermo, T, phase)) isa FT
            Test.@test Test.@inferred(MM2015.latent_heat(thermo, T, phase)) isa FT
            Test.@test Test.@inferred(MM2015.saturation(thermo, T, p, phase)) isa MM2015.PhaseSaturation{FT}
            Test.@test allocated(MM2015.saturation, thermo, T, p, phase) == 0
        end
        Test.@test Test.@inferred(MM2015.cp_m(thermo, q_t, q_l, q_i)) isa FT
        Test.@test Test.@inferred(MM2015.air_density(thermo, T, p, q_t, q_l, q_i)) isa FT
        Test.@test allocated(MM2015.air_density, thermo, T, p, q_t, q_l, q_i) == 0
        Test.@test allocated(MM2015.cp_m, thermo, q_t, q_l, q_i) == 0
    end
end
